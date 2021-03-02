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
#' @param distribution The distribution of the data. Used to determine in 
#'     which way to compute zscores. Either 'nb' or 'gaussian'.
#' @param peerResiduals If TRUE, PEER residuals are used to compute Z scores
#' @param ... Currently not used.
#' @return An OutriderDataSet containing the Z score matrix "zScore" as an 
#'     asasy.
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
#' ods <- computeEffectSizes(ods, effect_types=c("fold_change", "zscores"))
#' 
#' zScore(ods)[1:10,1:10]
#' assay(ods, "l2fc")[1:10,1:10]
#' 
#' @exportMethod computeZscores
setGeneric("computeZscores", 
        function(ods, ...) standardGeneric("computeZscores"))

#' #' @rdname computeZscores
#' #' @export
#' setMethod("computeZscores", "OutriderDataSet", 
#'     function(ods, peerResiduals=FALSE, ...){ 
#'         ZscoreMatrix(ods, peerResiduals=peerResiduals) })

#' @rdname computeZscores
#' @export
ZscoreMatrix <- function(ods, distribution=c("nb", "gaussian"), 
                        peerResiduals=FALSE, ...){ 
    
    if(length(normalizationFactors(ods)) == 0){
        stop("Please fit the autoencoder first before computing Z scores.")
    }
    
    distribution <- tolower(distribution)
    distribution <- match.arg(distribution)
    if(distribution == "nb"){
        residuals <- log2fc(ods)
    } else {
        residuals <- preprocessed(ods) - normalizationFactors(ods)
    }
    
    # Use residuals from PEER if present
    if(isTRUE(peerResiduals)){
        if(!"PEER_model" %in% names(metadata(ods)) && 
                !"residuals" %in% names(metadata(ods)[['PEER_model']])){
            stop("Please fit the data with 'peer' first.")
        }
        residuals <- metadata(ods)[['PEER_model']][['residuals']]
    }
    
    # default Zscore calculation
    Zscore <- (residuals - rowMeans(residuals, na.rm=TRUE)) / rowSds(residuals, 
                                                                    na.rm=TRUE)
    zScore(ods) <- Zscore
    validObject(ods)
    return(ods)
}

#' @rdname computeZscores
#' @export
setMethod("computeZscores", "Outrider2DataSet", ZscoreMatrix)
    
#' @param effect_types The types of effects to compute. Possible options are
#'     "fold_change" for log2 fold-changes, "zscores" for z scores, and 
#'     "delta" for the difference between observed and expected values.
#' @return An OutriderDataSet containing the requested effects as assays.
#'     Z score matrix in "zScore", the log2 fold changes in "l2fc", and/or 
#'     delta values in "delta".
#' 
#' @rdname computeZscores
#' @export
computeEffectSizes <- function(ods, distribution=c("nb", "gaussian"),
                        effect_types=c("fold_change", "zscores", "delta"), 
                        peerResiduals=FALSE){
    effect_types <- match.arg(effect_types, several.ok=TRUE)    
    distribution <- tolower(distribution)
    distribution <- match.arg(distribution)
    
    if("fold_change" %in% effect_types){
        if(any(preprocessed(ods) < 0, na.rm=TRUE)){
            warning("fold change calculation not meaninful for negative ",
                    "values, therefore computation is skipped. Use zscores ",
                    "and/or delta values instead.")
        } else{
            pc <- ifelse(distribution == "nb", 1, 0.01)
            assay(ods, "l2fc", withDimnames=FALSE) <- 
                log2fc(ods, pseudocount=pc)
        }
    }
    
    if("zscores" %in% effect_types){
        ods <- computeZscores(ods, distribution=distribution, 
                                    peerResiduals=peerResiduals)
    }
    
    if("delta" %in% effect_types){
        delta <- preprocessed(ods) - expected(ods)
        assay(ods, "delta", withDimnames=FALSE) <- delta
    }
    
    return(ods)
}


log2fc <- function(object, pseudocount=1){
    log2(preprocessed(object) + pseudocount) - 
        log2(normalizationFactors(object) + pseudocount)
}
