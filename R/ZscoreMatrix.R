#' 
#' Z score matrix
#' 
#' computes the z scores for every log2 fold-change compared
#' to the mean(normalizedCounts)
#' 
#' @param ods OutriderDataSet
#' @param ... further arguments passed to \code{ZscoreMatrix}
#' @return matrix of z scores and l2fc
#' 
#' @docType methods
#' @name computeZscores
#' @rdname computeZscores
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods <- estimateSizeFactors(ods)
#' ods <- computeZscores(ods)
#' 
#' assays(ods)[['zScore']][1:10,1:10]
#' assays(ods)[["l2fc"]][1:10,1:10]
#' 
#' @export
setGeneric("computeZscores", 
        function(ods, ...) standardGeneric("computeZscores"))

#' @rdname computeZscores
#' @export
setMethod("computeZscores", "OutriderDataSet", function(ods, ...){
    ZscoreMatrix(ods, ...)
})

ZscoreMatrix <- function(ods, normalized=TRUE, median=FALSE){
    log2fc <- log2fc(ods)
    Zscore <- (log2fc - rowMeans(log2fc)) / rowSds(log2fc)
    assays(ods)[["l2fc"]] <- log2fc
    assays(ods)[["zScore"]] <- Zscore
    validObject(ods)
    return(ods)
}


log2fc <- function(object, normalized = TRUE){
    l2fc <- log2(counts(object, normalized = normalized, offset = 1)) - 
        log2(rowMeans(counts(object, normalized = normalized, offset = 1)))
    return(l2fc)
}
    