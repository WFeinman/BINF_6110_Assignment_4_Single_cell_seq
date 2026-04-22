# Adjusted from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

#Presto package used for faster implementation of the Wilcoxon Rank Sum Test in seurat's FindMarkers function. Automatically used by seurat if installed without need of another library command.
#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')

#Package loading
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)

library(clusterProfiler)
library(org.Mm.eg.db)


library(DESeq2)
library(DOSE)

library(enrichplot)
library("dplyr")
library("GenomicFeatures")


library(tidyverse)


# Load the seurat dataset provided for this assignment
MouseSrt <- LoadSeuratRds("seurat_ass4.rds")

MouseSrt

# The [[ operator can add columns to object metadata. This is used to stash QC stats.
# Note: "pattern" is used to set a regex pattern to match features against. The "^mt" pattern works for mouse genomes, as per: https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/Seurat_QC_to_Clustering/
MouseSrt[["percent.mt"]] <- PercentageFeatureSet(MouseSrt, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(MouseSrt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Based on those plots, subsetting to mitochondrial percentage under 16% and # of genes over 600
MouseSrt <- subset(MouseSrt, subset = nFeature_RNA > 600 & percent.mt < 16)
VlnPlot(MouseSrt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MouseSrt


# Basic log normalization of the data (Could also instead use SCTransform()) and scaling
MouseSrt <- NormalizeData(MouseSrt, normalization.method = "LogNormalize")
MouseSrt <- FindVariableFeatures(MouseSrt, selection.method = "vst", nfeatures = 2000)

# By default Seurat only scales variable features. Here, we're instead scaling all features (better downstream visualization)
# The scaling phase is also where we would regress out unwanted sources of variation, e.g. cell cycle stage
#Note: Scaledata takes a while (~20-30 minutes), and a lot of memory (>30Gb).
all.genes <- rownames(MouseSrt)
MouseSrt <- ScaleData(MouseSrt, features = all.genes)

# We run a PCA to produce principal components that can be used to cluster our cells
MouseSrt <- RunPCA(MouseSrt, features = VariableFeatures(object = MouseSrt))


# Determine the ‘dimensionality’ of the dataset
ElbowPlot(MouseSrt)

# Clustering. Based on elbow plot, selecting 20 dimensions
MouseSrt <- FindNeighbors(MouseSrt, dims = 1:20)
MouseSrt <- FindClusters(MouseSrt, resolution = 0.5)

# UMAP v2 with clusters displayed
MouseSrt <- RunUMAP(MouseSrt, dims = 1:20)
DimPlot(MouseSrt, reduction = "umap", label = TRUE)

#Initial results seeem slighlty oversensitve to cluster mixing. Retrying with lower resolution
MouseSrt <- FindClusters(MouseSrt, resolution = 0.4)

# UMAP v2 with clusters displayed
MouseSrt <- RunUMAP(MouseSrt, dims = 1:20)
DimPlot(MouseSrt, reduction = "umap", label = TRUE)


# Cluster annotation
# find all markers of cluster 1
cluster1.markers <- FindMarkers(MouseSrt, ident.1 = 1)
head(cluster1.markers, n = 5)

# Top reported markers for cluster 1: "C1qa", "C1qb", "C1qc", "Ms4a7", and "Ctss"
VlnPlot(MouseSrt, features = c("C1qa", "C1qb", "C1qc", "Ms4a7", "Ctss"))

#FeaturePlot() is useful to display the features on our UMAP
FeaturePlot(MouseSrt, features = c("C1qa", "C1qb", "C1qc", "Ms4a7", "Ctss"))

#Since C1qa, C1qb, C1qc have  high overlap, mapping alternate features:
head(cluster1.markers, n = 10)
#Next 5 top-markers: "Cx3cr1" ,"Fcer1g", "Aif1", "Tyrobp", "Csf1r"
FeaturePlot(MouseSrt, features = c("Cx3cr1", "Fcer1g", "Aif1", "Tyrobp", "Csf1r"))


#C1qa, C1qb, and C1qc all code parts of the C1q subcomplex, which forms the C1 complex with C1r and C1s.
FeaturePlot(MouseSrt, features = c("C1ra", "C1rb", "C1s1", "C1s2"))



ifnb <- RunUMAP(MouseSrt, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("disease__ontology_label"))
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("time"))
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("organ_custom"))
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("mouse_id"))


#Take subset of data to focus analysis on cluster 1.
ifnb <- MouseSrt


# The actual aggregation
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("disease__ontology_label", "time", "organ_custom", "seurat_clusters"))


# DE test after pseudobulking our data
Idents(pseudo_ifnb) <- "seurat_clusters"

bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
                            ident.1 = "1",
                            ident.2 = "2",
                            test.use = "DESeq2")
head(bulk.mono.de, n = 15)



###Functional Annotation Section

# Map to Entrez IDs
res_df <- as.data.frame(bulk.mono.de)

ORF_ids <- rownames(res_df)


gene_map <- bitr(ORF_ids, 
                 fromType = "SYMBOL", 
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db)


# Add the mapping to results
res_df$SYMBOL <- sub("\\..*", "", rownames(res_df))
res_df <- merge(res_df, gene_map, by = "SYMBOL", all.x = TRUE)


# Define significant genes (Here we chose p < 0.05 and foldchange > 2, but other values could be justified)
sig_genes <- res_df %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()


# Define our background list of genes to compare to
# Rembember that ORA needs an "interesting gene" set and a background set.
all_genes <- res_df %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()



# Here we'll do a GO analysis for only Biological Process 

ego_bp <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff  = 0.2,
                   readable = TRUE)
head(as.data.frame(ego_bp))





# Now we'll do a KEGG analysis
kegg_enrich <- enrichKEGG(gene = sig_genes,
                          organism = 'mmu',
                          pvalueCutoff = 0.05,
                          qvalueCutoff  = 0.2)

head(as.data.frame(kegg_enrich))




# Dot plots
dotplot(ego_bp, showCategory = 20, title = "GO Biological Process")

dotplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment")


# Bar plot
barplot(ego_bp, showCategory = 15, title = "GO Biological Process")

barplot(kegg_enrich, showCategory = 15, title = "KEGG Biological Process")

# Enrichment map (Displays linked GO terms)
emapplot(pairwise_termsim(ego_bp), showCategory = 30)

# Notice that none of these plots split our genes by upregulated/downregulated?
# We would need to split those out ourselves as our gene set of interest
# Here's an example for upregulation:

upregulated_genes <- res_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

ego_bp_up <- enrichGO(gene = upregulated_genes,
                      universe = all_genes,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = FALSE)


#Plots of upregulated genes:

dotplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes from PC1 to PC2")

barplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes from PC1 to PC2")


#Kegg plot of upregulated genes


kegg_up <- enrichKEGG(gene = upregulated_genes,
                      organism = 'mmu',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)


barplot(kegg_up, showCategory = 15, title = "KEGG - Upregulated Genes from PC1 to PC2")

dotplot(kegg_up, showCategory = 15, title = "KEGG - Upregulated Genes from PC1 to PC2")


#Downregulate gene plots
downregulated_genes <- res_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC < -1) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

ego_bp_down <- enrichGO(gene = downregulated_genes,
                        universe = all_genes,
                        OrgDb = org.Mm.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2,
                        readable = FALSE)


dotplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes from PC1 to PC2")

barplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes from PC1 to PC2")


# as above, but for KEGG downregulation


kegg_down <- enrichKEGG(gene = downregulated_genes,
                        organism = 'mmu',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)


barplot(kegg_down, showCategory = 15, title = "KEGG - Downregulated Genes from PC1 to PC2")

dotplot(kegg_down, showCategory = 15, title = "KEGG - Downregulated Genes from PC1 to PC2")


#Testing Leukocyte marker (Ptprc) and T cell -markers Cd3d-g:
FeaturePlot(MouseSrt, features = c("Ptprc",  "Cd3d", "Cd3e", "Cd3g"))

#Leukocyte yest, t cell no. 
#Lymphocyte test? Itgam. pDCs test: Bst2. cDCs test? Itgax. Macrophage test? Adgre1. Monocyte test? Cd14.
FeaturePlot(MouseSrt, features = c("Itgam", "Bst2", "Itgax", "Adgre1", "Cd14"))


#Macrophage test. M1: Arg1, Cd38. M2: Egr2
FeaturePlot(MouseSrt, features = c("Arg1", "Cd38", "Egr2"))
