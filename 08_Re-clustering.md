# Re-clustering scRNA-seq data

## Overview
- Subset data to extract cluster or clusters for re-clustering
- Recalculate Variable Genes
- Scale
- Run PCA
- Re-cluster ( FindNeighbors + FindClusters)

## Codes
```
# df is the original seurat object that has gone through clustering, I'm interested to clusters 0-5 and 9, hypothetical neuron cells, for re-clustering
df_sub <- subset(df, idents = c(0, 1, 2, 3, 4, 5, 9))

# find out top 2000 variable genes
df_sub <- FindVariableFeatures(df_sub, selection.method = "vst", nfeatures = 2000)

# ScaleData by default will scale only 2000 variable genes, if using parameter features = all.gene, all genes will be scaled for heatmap visualization purpose, but still only 2000 variable genes will be used in the following RunPCA. 
all.genes <- rownames(df_sub)
df_sub <- ScaleData(df_sub, features = all.genes, verbose = FALSE)

df_sub <- RunPCA(df_sub, verbos = FALSE)
ElbowPlot(df_sub) 

df_sub <- RunPCA(df_sub, npcs = 13, verbose = FALSE) 

df_sub <- FindNeighbors(df_sub, reduction = "pca", dims = 1:13) 
df_sub <- FindClusters(df_sub, resolution = 0.5)

df_sub <- RunUMAP(df_sub, reduction = "pca", dims = 1:13) 
plot1 <- DimPlot(df_sub, reduction = "umap", label = TRUE, pt.size = 1)
pdf("control1.cluster.umap.recluster.pdf")
plot1 + coord_fixed()
dev.off()
```
