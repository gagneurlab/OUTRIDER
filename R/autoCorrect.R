#' 
#' autoCorrect
#' 
#' This is the wrapper function for the Python based autoCorrect method. 
#' It calls the autoencoder correction factor estimation.
#'
#' @param ods an OutriderDataSet object
#' @param predict run autoCorrect in predict mode, which means that the 
#'             prediction is based on a given model. Need a prefitted model.
#' @param save if TRUE, the fitted model will be saved to the given location.
#' @param epochs The number of epochs used for the autoCorrection training.
#' @param modelName Name of the model to read/write
#' @param modelDirectory The directory where the model is located 
#'             or should be stored
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
#' @export
autoCorrect <- function(ods, save=FALSE, predict=FALSE, epochs=250, 
                    modelName=NULL, modelDirectory=NULL){
    if(is.null(sizeFactors(ods))){
        stop(paste("Please calculate the size factors before calling", 
                "the autoCorrect function"))
    }
    if(is.null(modelName)){
        modelName <- "autoCorrectModel"
    }
    if(is.null(modelDirectory)){
        modelDirectory <- file.path(tempdir(), "models")
    }
    # try loading the autoCorrect module
    if(!checkAutoCorrectExists()){
        warning(paste("Can not load the autoCorrect python module.", 
                "Please make sure that it is within the python path.", 
                "No correction was done!"))
        return(ods)
    }
    if(isScalarNumeric(options("OUTRIDER.epochs"))){
        epochs <- options("OUTRIDER.epochs")
    }
    epochs <- as.integer(epochs)
    
    # get needed data
    k <- counts(ods, normalized=FALSE)
    sf <- sizeFactors(ods)
    
    # transform data to be used by autoCorrect
    kt <- t(as.matrix(k))
    sfm <- as.matrix(sf)
    
    # correctionFactors is a matrix of the same dimension as k
    autoCorrectObj <- import("autoCorrection")
    corrected <- autoCorrectObj$correctors$AECorrector(modelName,
            modelDirectory, save_model=save, epochs=epochs)$correct(
                    kt, sfm, only_predict=predict)
    correctionFactors <- t(corrected)
    stopifnot(identical(dim(k), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    validObject(ods)
    return(ods)
}

#'
#' Checks if the autoCorrection python module exists
#' 
#' @noRd
checkAutoCorrectExists <- function(){
    try({autoCorrectObj <- import("autoCorrection")})
    exists("autoCorrectObj")
}

