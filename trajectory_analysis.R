#################First trajectory: Control - 0h - 12h - 24h#####################
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

gam <- qs_read('mgcGam.qs2')
res <- associationTest(gam)
sigRes <- subset(res, pvalue < 1e-2)
sigRes$df <- c()
sigRes <- addRanks(sigRes, signs = c(-1, 1, -1))
write.csv(sigRes, 'Gene associated with pseudotime - Ctrl-0h-12h-24h.csv')

a <- subset(res, pvalue < 1e-4 & meanLogFC > 1 & waldStat > 100)
sigGenes <- rownames(a)


featureWes(miniSeurat, 'FTH1', idClass = 'orig.ident')
DimPlot(miniSeurat, group.by='orig.ident')

#################Second trajectory: Control - 24h - 12h - 0h#####################
miniSeurat <- qs_read('MGCSeurat005.qs2')

sce <- as.SingleCellExperiment(miniSeurat)
sce <- slingshot(sce,
                 clusterLabels = 'orig.ident',
                 reducedDim = 'UMAP',
                 start.clus = 'Control',
                 end.clus = '0h',
                 dist.method = 'simple')

miniSeurat <- addLineages(miniSeurat, sce)
mat <- extractLineageMat(sce, 'Lineage1')
df <- points2Seg(mat)
p <- singleLineagePlot(miniSeurat, sce, 'Lineage1', 'orig.ident')

p1 <- featureWes(miniSeurat, 'Lineage1', idClass='orig.ident',
                 labelSize=4) + labs(color='Pseudotime') +
    geom_segment(data=df, aes(x=x, y=y, xend=xEnd, yend=yEnd),
                 arrow = arrow(length = unit(0.1, "cm")))

v <- subset(miniSeurat, features=VariableFeatures(miniSeurat)[1:100])
