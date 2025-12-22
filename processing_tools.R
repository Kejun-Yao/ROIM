basicDimRed <- function(seuratObj){
    DefaultAssay(seuratObj) <- 'RNA'
    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- FindVariableFeatures(seuratObj)
    seuratObj <- ScaleData(seuratObj)
    seuratObj <- RunPCA(seuratObj)
    
    DefaultAssay(seuratObj) <- 'ATAC'
    seuratObj <- FindTopFeatures(seuratObj, min.cutoff=5)
    seuratObj <- RunTFIDF(seuratObj)
    seuratObj <- RunSVD(seuratObj)
    
    DefaultAssay(seuratObj) <- 'RNA'
    return(seuratObj)
}

jointIntegration <- function(seuratObj, 
                             dimsList = list(1:50, 2:40),
                             seed = 1){
    seuratObj <- with_seed(seed, 
                           RunHarmony(seuratObj, 
                                      group.by.vars = 'orig.ident', 
                                      reduction.use = 'pca', 
                                      reduction.save = 'harmony_rna', 
                                      assay.use = 'RNA',
                                      project.dim = FALSE))
    seuratObj <- with_seed(seed,
                           RunHarmony(seuratObj, 
                                      group.by.vars = 'orig.ident', 
                                      reduction.use = 'lsi', 
                                      reduction.save = 'harmony_atac', 
                                      assay.use = 'ATAC',
                                      project.dim = FALSE))
    seuratObj <- FindMultiModalNeighbors(seuratObj, 
                                         reduction.list = list('harmony_rna', 
                                                               'harmony_atac'), 
                                         dims.list = dimsList,
                                         verbose = TRUE)
    seuratObj <- RunUMAP(seuratObj, nn.name = "weighted.nn", verbose=TRUE)
    DefaultAssay(seuratObj) <- 'RNA'
    return(seuratObj)
}