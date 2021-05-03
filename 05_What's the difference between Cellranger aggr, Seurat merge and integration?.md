# What's the difference between Cellranger aggr, Seurat merge and integration?
https://www.reddit.com/r/bioinformatics/comments/n3ltnf/cellranger_aggr_vs_seurat_merge_and_integrate/

## Cellranger aggr:
- input is raw count matrix
- Chemistry Batch Correction (different versions of kit reagents)
- algorithm is based on mutual nearest neighbors (MNN)
- mapped (default): For each library type, subsample reads from higher-depth GEM wells until they all have, on average, an equal number of reads per cell that are confidently mapped to the transcriptome (Gene Expression) or assigned to known features (Feature Barcode Technology). Can set to none. 

## Seurat merge
- merged raw count matrix. merges the raw count matrices of two Seurat objects and creates a new Seurat object with the resulting combined raw count matrix. 

## Seurat SCTransform
- can regress out batch effects etc

## Seurat integration
- can be used both to correct for technical differences between datasets (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions.

# What's the difference between batch effects and data integration
- batch effects between samples or cells in the same experiment is the classical scenario known as batch correction from bulk RNA-seq. (same experiment, different batches)
- The integration of data from multiple experiments, which we call data integration. (different experiment/conditions)
- batch correction method: ComBat
- data integration method: CCA, MNN, Scanorama, RISC, scGen, LIGER, BBKNN, Harmony. 
- While data integration methods can also be applied to simple batch correction problems.

# What's the difference between functions of ScaleData and SCTransform
In Seurat V2, the standard pipeline is: NormalizeData(), FindVariableFeatures() and ScaleData(). Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA. 
```
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```
Seurat V3 still has the ScaleData function, but new function SCTransform has replaced functions of NormalizeData(), FindVariableFeatures(), ScaleData(). During normalization, we can also remove confounding sources of variation.

```
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
```




