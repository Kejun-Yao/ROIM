library(henna)

featPatchwork <- function(seuratObj, features){
    plots <- lapply(features, function(x) 
        featureWes(seuratObj, feature=x, idClass='orig.ident', repel=TRUE))
    p <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
    return(p)
}

miniSeurat <- qs_read('MGCSeurat005MCA.qs2')
table1 <- read.csv('Results/Genes associated with pseudotime - Ctrl-0h-12h-24h.csv')
table2 <- read.csv('Results/Genes associated with pseudotime - Ctrl-24h-12h-0h.csv')

genes <- table1[seq(30), 1]
p <- genesDimPlot(miniSeurat, table1[seq(30), 1], groupBy='orig.ident') + 
    ggtitle('Centers of mass - top pseudotime-linked genes')
devPlot(p)

genes <- table1[seq(30), 1]
m <- scExpMat(miniSeurat, genes=genes)
m <- data.frame(cmdscale(dist(m)))
p <- densityPlot(m, 'MDS plot - Top pseudotime-linked genes', drawNN=FALSE)
devPlot(p)

p <- featPatchwork(miniSeurat, c('HSPB1', 'HSPA8', 'HSPH1', 'CRYAB'))
devPlot(p)

p <- featPatchwork(miniSeurat, c('FTL', 'FTH1', 'NEAT1', 'HSP90AB1'))
devPlot(p)

p <- featPatchwork(miniSeurat, c('MT-ATP6', 'MT-CYB', 'TF', 'NAV2'))
devPlot(p)

genes <- table2[seq(30), 1]
p <- genesDimPlot(miniSeurat, table1[seq(30), 1], groupBy='orig.ident') + 
    ggtitle('Centers of mass - top pseudotime-linked genes')
devPlot(p)

genes <- table2[seq(30), 1]
m <- scExpMat(miniSeurat, genes=genes)
m <- data.frame(cmdscale(dist(m)))
p <- densityPlot(m, 'MDS plot - Top pseudotime-linked genes', drawNN=FALSE)
devPlot(p)


