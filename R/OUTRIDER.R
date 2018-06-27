#' 
#' OUTRIDER - Full analysis pipeline
#' 
#' The OUTRIDER function runs the default OUTRIDER pipline. 
#' Combinig the fit, the computation of zScores and pValues.
#' All computed values are returned as a OutriderDataSet object.
#' 
#' @param object main OutriderDataSet object, which contains all the data 
#' @param autoCorrect if TRUE the raw read counts will be controled 
#'             for confounders by the autoencoder 
#' @param nbModelFile The file where the negative binomial model should be saved
#' @return OutriderDataSet with all the computed values. The values can be 
#'             accessed by: \code{assays(ods)[['value']]}
#'
#' @examples
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods <- OUTRIDER(ods)
#' 
#' assays(ods)[['pValue']][1:10,1:10]
#' res <- results(ods, all=TRUE)
#' plotVolcano(ods, 1)
#' 
#' @export
OUTRIDER <- function(object, autoCorrect=TRUE, nbModelFile=NULL){
    message(paste0(date(), ": SizeFactor estimation ..."))
    object <- estimateSizeFactors(object)
    if(autoCorrect == TRUE){
        message(paste0(date(), ": Running auto correct ..."))
        object <- autoCorrect(object, q=20)
    }
    if(is.null(nbModelFile) || !file.exists(nbModelFile)){
        message(paste0(date(), ": Fitting the data ..."))
        object <- fit(object, modelFile=nbModelFile)
        nbModelFile=NULL
    } else {
        message("Using existing NB parameters from model: ", nbModelFile)
    }
    message(paste0(date(), ": P-value calculation ..."))
    object <- computePvalues(object, modelFile=nbModelFile)
    message(paste0(date(), ": Zscore calculation ..."))
    object <- computeZscores(object)
    return(object)
}

