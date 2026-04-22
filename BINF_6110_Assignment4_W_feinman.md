---
noteID: 2c7484a1-3e13-4c86-9553-1880a078cd0f
---
# Single-Cell Transcript Mapping of Mouse Epithelial C1q Distribution during Influenza A infection.
## Introduction

Single cell sequencing is a technology with remarkable potential utility in a wide variety of biological contexts. Essentially, a slurry of cells are isolated from one or more sample tissues, isolated through a variety of capture mechanisms, suspended in a nano-droplet with a unique identifier attached. By cataloging the transcriptome of identifiable individual cells and mapping them, it is possible to not only distinguish tissue types by RNA output, but find specializations and subdivisions within tissue types, discern inter-cellular expression patterns, and gain a much better idea of exactly how specific cells respond to a given stimulus in a tissue context, rather than taking average estimations or extrapolating from an artificial monoculture. (Jovic et al., 2022)

A recent paper used this technology to better characterize pathology and immune response in mice infected with Influenza A. Their study provides a large dataset describing single cell expression changes in nasal epithelial tissue samples over the course of a two week infection. While their study touched on much, ranging from neutrophil production to a rare subset of epithelial floor cells, much more can potentially be gleaned from the dataset. (Kazer et al., 2024)

This code is exploratory in nature, taking a look at a discrete subset of the cells identified through Principal Component Analysis (PCA), mapping its RNA expression patterns, and attempting to draw wider conclusions about its role in the transcriptome through functional annotation.


## Methods

A Seurat object compiling the transcriptome data of the study was provided as part of this course: [https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds](https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds)

If reconstruction of the Seurat object is desired, the accession numbersfor the original study may be accessed on the key findings table here: https://pmc.ncbi.nlm.nih.gov/articles/PMC11324402/#SM1

Additional code used in the original study, such as metadata for the Seurat object, may be found here: https://doi.org/10.5281/zenodo.11104940

This paper primarily focused on use of the Seurat R software package to explore transcript data, setting initial parameters using an exploratory violin plot to filter outliers, then an elbow plot to determine optimal dimensionality for clustering. Once inputs were determined determined, PCA was run with varying resolution until a compromise was found between cluster fidelity and group distinctiveness. 

A distinct cluster with clearly identifiable boundaries and transition points to other clusters was selected. In this case, that was Cluster 1, transitioning to Cluster 2. Exploratory UMAP was performed to try and identify which genes were responsible for that transition, and any underlying explanation for that transition in the Seurat metadata.

Afterwards, the clusterProfiler package was used to characterize the transition from a biological perspective, using functional annotation to identify what processes and gene groups were involved in the transition.


## Results
### UMAP with Clusters (Resolution 0.5)
![[Pasted image 20260421131841.png]]

### UMAP with Clusters (Resolution 0.4)
![[Pasted image 20260421131509.png]]
**Figures 1 and 2:** UMAP clusters of cell populations by PCA. Figure 1 is included to show there is some degree of variation within each of clusters, though displays some over-mixing. Figure 2 shows a compromise preserving cluster identity while minimizing overmixing.

### Feature plots of top 10 genes by p_value for Cluster 1
![[Pasted image 20260421145504.png]]
![[Pasted image 20260421145516.png]]
**Figures 3 and 4:** Feature plots of 10 most significant markers of cluster 1, overlaid on UMAP. Intensity of purple indicates degree of expression. 

C1qa, C1qb, C1qc, and Ms4a7 were the most specific predictors for Cluster 1. In the classical C1 pathway, C1q combined with C1s and C1r, forming the C1 complex. This binds to antigen-antibody complexes and activates the complement system of the immune response. (Venkatraman Girija et al., 2013)

An alternative mechanism may be active, however. Without C1s and C1r, C1q will enhance FcγR and CR1 mediated phagocytosis in cases of low antibody coverage on the target. They also play a role in apoptic cell phagocytosis and clearance. (Thielens et al., 2017)

Ms4a7 is from the 4A gene membrane spanning family, components of varying signaling pathways across different tissue types. Its human ortholog MS4A7 contributes to inflammatory signaling. (Liang & Tedder, 2001; Zhou et al., 2025)

Ctss codes for Cathepsin S, a cysteine protease secreted during inflammation, and interacts with proinflammatory cytokines. One possibility for the overlap between C1qa-c, Ms4a7m and Ctss may be co-activation of immune response and inflammatory pathways.  (Ainscough et al., 2017)

A follow-up mapping of other C1 components was made to determine if the C1 pathway was being upregulated.
### Feature plots of other C1 complex components: 
![[Pasted image 20260421161540.png]]
**Figure 5:** Feature plot of C1ra, C1s1, C1s2. C1rb plotting was attempted, but no match found.

Significantly lower expression levels of other components of the C1 complex were displayed compared to C1qa-c. Additionally, what components were displayed were isolated to different clusters, indicating the C1 complex mechanism is not actively transcribed in Cluster 1, and the C1q sub-complex was being transcribed independently.

To clarify this shift between PC1 and PC2 in C1qa-c and Ms4a7, mapping of alternative explanations via sample metadata were tested.

### UMAP Plots of Sample Types

![[Pasted image 20260421180922.png]]
**Figure 6:** Annotated UMAP of influenza infection status.

![[Pasted image 20260421181941.png]]
**Figure 7:** Annotated UMAP of infection timeline. Label indicates days-post-infection of Influenza A virus.

![[Pasted image 20260421182047.png]]
**Figure 8:** UMAP of tissue origin site. LNG = Lateral Nasal Gland. OM = Olfactory Mucosa. RM = Respiratory Mucosa.

Overall, Seurat metadata was inconclusive. Some trends exist (such as the greater prevalence of RM in cluster 2), but not enough to explain the transition in an of itself.
### Enrichment Map: Transcript Activity by Linked GO terms
![[Pasted image 20260421232636.png]]
**Figure 9:** Enrichment map by GO-term linkage. Shows functional relations of altered gene expression between PC1 and PC2 in a wider biological context, grouping related terms together.

![[Pasted image 20260421233217.png]]
![[Pasted image 20260421233314.png]]
**Figures 10 and 11:** Clarify upregulated activities in PC2 relative to PC1 in a gene ontology (Figure 10) and pathway context (Figure 11).

Overall upregulation trends towards PC2: Increased leukocyte migration and chemotaxis. Cytokine-receptor interaction. Viral protein interaction with cytokines/receptors. MAPK/ERK signaling pathway activity. Lysosome biogenesis.

![[Pasted image 20260421234920.png]]
![[Pasted image 20260421234953.png]]
**Figures 12 and 13**: Clarify downregulated activities in PC2 relative to PC1 in a gene ontology (Figure 12) and pathway context (Figure 13).

Overall downregulation trends towards PC2: Lymphocyte activation regulation. Cell-adhesion, particularly in leukocytes. Leukocyte proliferation. Hallmarks of herpes virus infection. T cell proliferation.

### UMAP Feature Plots of Lymphocyte Markers:
![[Pasted image 20260422003912.png]]
**Figure 14:**  UMAP of Lymphocyte signaling markers. Itgam is a broad Lymphocyte marker. Bst2 marks for Plasmacytoid Dendritic Cells. Itgax marks for Conventional Dendritic Cells. Adgre1 marks for Macrophages. Cd14 marks for Monocytes. (NCBI, 2026)

![[Pasted image 20260422012352.png]]
**Figure 15:** Other Macrophage M1 and M2 distinguishing markers. Cd38 labels M1, Arg1 and Egr2 label M2. (Jablonski et al., 2015)

While Adgre1 indicated a strong macrophage signal, subsequent macrophage-specific markers have only offered partial confirmation.

## Discussion

Overall, as is somewhat characteristic of single cell sequencing, this analysis posed more questions than it answered. It's a credit to the technique that multiple rounds of investigation into the data produced can continue to provide further investigation avenues.

On a general note, as with the original study, the transition between human and mice orthologs of a given gene means some qualifiers must be noted. While both humans and mice have been extensively studied, those fields of study do not perfectly overlap. Gene functions characterized in a human may not necessarily be documented or be applicable in mouse models, or vice versa. This is of particular note when trying to extrapolate gene functions from mice with only a tangential human ortholog to draw upon, especially in highly variable, diverse function protein families. (Kazer et al., 2024)

The primary identifiers of cluster 1 were C1qa, C1qb, C1qc, and Ms4a7, displayed constitutively throughout the cluster and forming a transition series into PC2. C1qa-c is an integral part of the PC1 complex, bonding alongside the C1r and C1s proteins to form the classical complement pathway, attaching to antibody-antigen complexes and recruiting immune cells to the target site. The collagen-like tail of the C1q subcomplex forms a vital attachement site for this purpose. (Sellar et al., 1992; Thielens et al., 2017; Venkatraman Girija et al., 2013)

However, the other components of the classical complement pathway were notably absent in the region. C1s and C1r orthologs showed negligible expression in PC1. What is likely occurring, therefore, is an alternative pathway independent of the C1 complex. C1q is capable of enhancing phagocytosis in response to antibody-poor antigen targets and of flagging apoptotic cells. (Thielens et al., 2017)

Ms4a7 showed a promising degree of overlap with C1q, indicate correlation with it's expression. However, the functional interactions of the mouse ortholog of the signaling protein have not been thoroughly investigated in research literature. The human ortholog MS4A7 has been demonstrated to play a role in inflammation-response, which would correlate with immune response from C1q and the wider circumstances of viral infection, but that is extrapolation across species. (Liang & Tedder, 2001; Zhou et al., 2025)

Likewise, functional annotation continued to pose questions. The MAPK signaling cascades are another potential avenue for signal correlation though like the MS4A family they serve a wide variety of functions. (Plotnikov et al., 2011)

Lymphocyte activity in the area would correlate with immune flagging from C1q, and some evidence to support that was shown in immune marker profiles. However, of the lymphocytes identified as displating activity in the area, identifying one cell type with perfect correlation has proven troublesome. The best correlator identified in that regard was Adgre1, a macrophage marker, but subsequent markers for M1 and M2 were unable to show as clear a correlation. (Jablonski et al., 2015; Li et al., 2026)

The functional annotation reporting on herpes related activity in PC1 relative to PC2 was unexpected. As their was already an Influenza A infection occurring in the sample, common viral immune responses may have been flagged as herpes, though it is also possible a latent herpes infection had not been explicitly tested for.

Overall, it would appear C1q subunits and Ms4a7 have a correlative role as part of inflammation response in mouse epithelial cells, though further mechanistic testing would need to be done to ascertain exactly what that relationship is.


### References

Ainscough, J. S., Macleod, T., McGonagle, D., Brakefield, R., Baron, J. M., Alase, A., Wittmann, M., & Stacey, M. (2017). Cathepsin S is the major activator of the psoriasis-associated proinflammatory cytokine IL-36γ. _Proceedings of the National Academy of Sciences_, _114_(13). [https://doi.org/10.1073/pnas.1620954114](https://doi.org/10.1073/pnas.1620954114)

Ferretti, E., Pistoia, V., & Corcione, A. (2014). Role of Fractalkine/CX3CL1 and Its Receptor in the Pathogenesis of Inflammatory and Malignant Diseases with Emphasis on B Cell Malignancies. _Mediators of Inflammation_, _2014_, 1–10. [https://doi.org/10.1155/2014/480941](https://doi.org/10.1155/2014/480941)

Jablonski, K. A., Amici, S. A., Webb, L. M., Ruiz-Rosado, J. D. D., Popovich, P. G., Partida-Sanchez, S., & Guerau-de-Arellano, M. (2015). Novel Markers to Delineate Murine M1 and M2 Macrophages. _PLOS ONE_, _10_(12), e0145342. [https://doi.org/10.1371/journal.pone.0145342](https://doi.org/10.1371/journal.pone.0145342)

Jovic, D., Liang, X., Zeng, H., Lin, L., Xu, F., & Luo, Y. (2022). Single‐cell RNA sequencing technologies and applications: A brief overview. _Clinical and Translational Medicine_, _12_(3), e694. [https://doi.org/10.1002/ctm2.694](https://doi.org/10.1002/ctm2.694)

Kazer, S. W., Match, C. M., Langan, E. M., Messou, M.-A., LaSalle, T. J., O’Leary, E., Marbourg, J., Naughton, K., Von Andrian, U. H., & Ordovas-Montanes, J. (2024). Primary nasal influenza infection rewires tissue-scale memory response dynamics. _Immunity_, _57_(8), 1955-1974.e8. [https://doi.org/10.1016/j.immuni.2024.06.005](https://doi.org/10.1016/j.immuni.2024.06.005)

Li, J., Lu, J., Zheng, C., Huang, X., Li, H., Mai, Q., Chen, S., Zhou, Z., Zhu, J., Yu, T., Xu, M., Tan, H., Zhang, C., Gao, Q., Liu, J., & Pan, C. (2026). Serotonin-licensed macrophages potentiate chemoresistance via inositol metabolic crosstalk in ovarian cancer. _Cell Metabolism_, _38_(2), 331-349.e10. [https://doi.org/10.1016/j.cmet.2025.11.011](https://doi.org/10.1016/j.cmet.2025.11.011)

Liang, Y., & Tedder, T. F. (2001). Identification of a CD20-, FcϵRIβ-, and HTm4-Related Gene Family: Sixteen New MS4A Family Members Expressed in Human and Mouse. _Genomics_, _72_(2), 119–127. [https://doi.org/10.1006/geno.2000.6472](https://doi.org/10.1006/geno.2000.6472)

Plotnikov, A., Zehorai, E., Procaccia, S., & Seger, R. (2011). The MAPK cascades: Signaling components, nuclear roles and mechanisms of nuclear translocation. _Biochimica et Biophysica Acta (BBA) - Molecular Cell Research_, _1813_(9), 1619–1633. [https://doi.org/10.1016/j.bbamcr.2010.12.012](https://doi.org/10.1016/j.bbamcr.2010.12.012)

Sellar, GrantC., Cockburn, D., & Reid, KennethB. M. (1992). Localization of the gene cluster encoding the A, B, and C chains of human C1q to 1p34.1?1p36.3. _Immunogenetics_, _35_(3). [https://doi.org/10.1007/BF00185116](https://doi.org/10.1007/BF00185116)

Thielens, N. M., Tedesco, F., Bohlson, S. S., Gaboriaud, C., & Tenner, A. J. (2017). C1q: A fresh look upon an old molecule. _Molecular Immunology_, _89_, 73–83. [https://doi.org/10.1016/j.molimm.2017.05.025](https://doi.org/10.1016/j.molimm.2017.05.025)

Venkatraman Girija, U., Gingras, A. R., Marshall, J. E., Panchal, R., Sheikh, Md. A., Harper, J. A. J., Gál, P., Schwaeble, W. J., Mitchell, D. A., Moody, P. C. E., & Wallis, R. (2013). Structural basis of the C1q/C1s interaction and its central role in assembly of the C1 complex of complement activation. _Proceedings of the National Academy of Sciences_, _110_(34), 13916–13920. [https://doi.org/10.1073/pnas.1311113110](https://doi.org/10.1073/pnas.1311113110)

Zhou, L., Lu, Y., Qiu, X., Chen, Z., Tang, Y., Meng, Z., Yan, C., Du, H., Li, S., & Lin, J. D. (2025). Lipid droplet efferocytosis attenuates proinflammatory signaling in macrophages via TREM2- and MS4A7-dependent mechanisms. _Cell Reports_, _44_(2), 115310. [https://doi.org/10.1016/j.celrep.2025.115310](https://doi.org/10.1016/j.celrep.2025.115310)

