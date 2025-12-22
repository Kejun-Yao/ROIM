jointReadSeurats <- function(folder, folderPrefix, ident, annotation){
    folderPath <- file.path(folderPrefix, folder, 'outs')
    counts <- Read10X(file.path(folderPath, 'filtered_feature_bc_matrix'))
    fragpath <- file.path(folderPath, 'atac_fragments.tsv.gz')
    seuratObj <- CreateSeuratObject(
        counts = counts$`Gene Expression`,
        project = ident,
        assay = "RNA"
    )
    seuratObj[["ATAC"]] <- CreateChromatinAssay(
        counts = counts$Peaks,
        sep = c(":", "-"),
        fragments = fragpath,
        annotation = annotation
    )
    qs_save(seuratObj, paste0(folder, 'RawSeurat.qs2'))
    return(seuratObj)
}