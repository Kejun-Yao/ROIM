miniSeurat <- qs_read('MGCSeurat005.qs2')

sce <- as.SingleCellExperiment(miniSeurat)
sce <- slingshot(sce,
                 clusterLabels = 'orig.ident',
                 reducedDim = 'UMAP',
                 start.clus = 'Control',
                 dist.method = 'simple')

miniSeurat <- addLineages(miniSeurat, sce)
mat <- extractLineageMat(sce, 'Lineage1')
df <- points2Seg(mat)
p <- singleLineagePlot(miniSeurat, sce, 'Lineage1', 'orig.ident')

p1 <- featureWes(miniSeurat, 'Lineage1', idClass='orig.ident',
                 labelSize=4) + labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))

pseudotime <- slingPseudotime(sce)
cellWeights <- slingCurveWeights(sce)

counts <- LayerData(miniSeurat, assay='RNA', layer='counts')

x <- Sys.time()
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 60

gam <- fitGAM(counts,
              pseudotime=pseudotime,
              cellWeights=cellWeights,
              parallel=TRUE,
              genes=VariableFeatures(miniSeurat),
              BPPARAM=BPPARAM)
qs_save(tsGAM, 'mgcGamVarGenes.qs2')

y <- Sys.time()
print(y - x)


