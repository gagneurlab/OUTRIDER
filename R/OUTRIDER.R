#' 
#' The OUTRIDER function runs the default OUTRIDER pipline. 
#' Combinig the fit, the computation of zScores and pValues.
#' All computed values are returned as a OutriderDataSet object.
#' 
#' The values can be accessed by: assays(ods)[['value']] 
#' 
#' @param object main OutriderDataSet object, which contains all the data 
#' @param autoControl option, to run autoControl
#'
#' @return OutriderDataSet with all the computed values
#'
#' @examples
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods <- OUTRIDER(ods)
#'   
#' assays(ods)[['pValue']][1:10,1:10]
#' res <- results(ods)
#' plotVolcano(ods, 1)
#'   
#' @export
OUTRIDER <- function(object, autoControl=TRUE, save=FALSE,
            acModelName=NULL, acModelDir=NULL, nbModelFile=NULL){
    message(paste0(date(), ": SizeFactor estimation ..."))
    object <- estimateSizeFactors(object)
    if(autoControl == TRUE){
        message(paste0(date(), ": Running auto correct ..."))
        predict=FALSE
        if(!is.null(acModelName) & !is.null(acModelDir) & save == FALSE){
            predict=TRUE
        }
        tryCatch({object <- autoCorrect(object, save=save, predict=predict, 
            modelName=acModelName, modelDirectory=acModelDir)},
                error=function(e) warning(e))
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