#' 
#' Z score computation
#' 
#' Computes the z scores for every count in the matrix.
#' The z score is defined in the log\eqn{_2}{[2]} space as follows:
#' \eqn{z_{ij} = \frac{l_{ij} - \mu_j^l}{
#'      \sigma_j^l}}{z_ij = (l[ij] - mu[j]^l)/sigma[ij]}
#' where \code{l} is the log\eqn{_2}{[2]} transformed normalized count and 
#' \eqn{\mu}{mu} and \eqn{\sigma}{sigma} the mean and standard deviation 
#' for gene \code{j} and sample \code{i}, respectively.
#' 
#' @param ods OutriderDataSet
#' @param ... Further arguments passed on to \code{ZscoreMatrix}.
#' @param peerResiduals If TRUE, PEER residuals are used to compute Z scores
#' @return An OutriderDataSet containing the Z score matrix "zScore" and
#'     the log2 fold changes "l2fc" as asasys.
#' 
#' @docType methods
#' @name computeZscores
#' @rdname computeZscores
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' 
#' ods <- controlForConfounders(ods, implementation="pca")
#' ods <- computeZscores(ods)
#' 
#' zScore(ods)[1:10,1:10]
#' assay(ods, "l2fc")[1:10,1:10]
#' 
#' @exportMethod computeZscores
setGeneric("computeZscores", 
        function(ods, ...) standardGeneric("computeZscores"))

#' @rdname computeZscores
#' @export
setMethod("computeZscores", "OutriderDataSet", 
    function(ods, peerResiduals=FALSE, ...){ 
        ZscoreMatrix(ods, peerResiduals=peerResiduals) })

ZscoreMatrix <- function(ods, peerResiduals){
    if(length(normalizationFactors(ods)) == 0){
        stop("Please fit the autoencoder first befor computing Z scores.")
    }
    
    # default Zscore calculation
    log2fc <- log2fc(ods)
    Zscore <- (log2fc - rowMeans(log2fc)) / rowSds(log2fc)
    
    # Use residuals from PEER if present
    if(isTRUE(peerResiduals)){
        if(!"PEER_model" %in% names(metadata(ods)) && 
                    !"residuals" %in% names(metadata(ods)[['PEER_model']])){
            stop("Please fit the data with 'peer' first.")
        }
        residuals <- metadata(ods)[['PEER_model']][['residuals']]
        Zscore <- (residuals - rowMeans(residuals)) / rowSds(residuals)
    }
    
    assay(ods, "l2fc", withDimnames=FALSE) <- log2fc
    zScore(ods) <- Zscore
    validObject(ods)
    return(ods)
}


log2fc <- function(object){
    log2(counts(object) + 1) - log2(normalizationFactors(object) + 1)
}

