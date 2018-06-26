#' 
#' Z score computation
#' 
#' Computes the z scores for every count in the matrix. 
#' The z score is defined in the log2 space as follows:
#' \ifelse{html}{
#'     \out{z<sub>ij</sub> = (l<sub>ij</sub> - mu<sub>j</sub><sup>l</sup>)/
#'             sigma<sub>j</sub><sup>l</sup>}}{
#'     \deqn{z_{ij} = \frac{l_{ij} - \mu_j^l}{\sigma_j^l}}},
#' where l is the log2 transformed normalized count and mu and sigma the 
#' mean and standard deviation for gene j, respectively.
#' 
#' @param ods OutriderDataSet
#' @param ... Further arguments passed to \code{ZscoreMatrix}
#' @return An OutriderDataSet containing the z score matrix ("zScore") and
#'     the log2 fold changes ("l2fc") as asasys.
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
#' @exportMethod computeZscores
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
    