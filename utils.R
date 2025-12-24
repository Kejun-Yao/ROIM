addRanks <- function(df, 
                     columns = colnames(df), 
                     signs = rep(1, ncol(df)),
                     tiesMethod = 'min',
                     useMeanRank = TRUE,
                     replaceByMeanRank = TRUE){
    rankColumns <- paste0(columns, 'Rank')
    for (i in seq_along(columns)){
        column <- columns[i]
        sign <- signs[i]
        rankColumn <- rankColumns[i]
        df <- df[order(df[[column]]), ]
        df[[rankColumn]] <- rank(sign * df[[column]], 
                                 ties.method=tiesMethod)
    }
    if (useMeanRank){
        df$rank <- rowMeans(df[, rankColumns])
        if(replaceByMeanRank)
            df[, rankColumns] <- c()
        
    }
    df$rank <- rank(df$rank, ties.method=tiesMethod)
    df <- df[order(df$rank), ]
    return(df)
}