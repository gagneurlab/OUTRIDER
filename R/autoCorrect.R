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
#' @param param_path The directory where the optimal parameters are stored
#' @param param_exp_name Name of the parameter file
#' @param verbose if TRUE further information about the training/predicting
#'             of the autoencoder is printed.
#' @param seed if this is an integer it will be used to set the seed in the 
#'             python runtime environment
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' ods <- autoCorrect(ods, seed=42)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
#' @export
autoCorrect <- function(ods, save=FALSE, predict=FALSE, epochs=250, 
                    param_path=NULL, param_exp_name=NULL, verbose=FALSE,
                    modelName=NULL, modelDirectory=NULL, seed=NULL){
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
    stopifnot(isScalarLogical(verbose))
    verbose <- as.integer(verbose)
    if(is.numeric(seed)){
        seed <- as.integer(seed)
    }
    
    # get needed data
    k <- counts(ods, normalized=FALSE)
    sf <- sizeFactors(ods)
    
    # transform data to be used by autoCorrect
    kt <- t(as.matrix(k))
    sfm <- as.matrix(sf)
    
    # correctionFactors is a matrix of the same dimension as k
    autoCorrectObj <- import("autoCorrection")
    corrected <- autoCorrectObj$correctors$AECorrector(epochs=epochs, modelName,
            modelDirectory, save_model=save, epochs=epochs, verbose=verbose,
            param_exp_name=param_exp_name, param_path=param_path, 
            seed=seed)$correct(
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
    if(!exists("autoCorrectObj")){
        return(FALSE)
    }
    
    # Check if autoCorrection version is correct
    aeVersion <- as.character(py_get_attr(autoCorrectObj, '__version__'))
    if(compareVersion(aeVersion, "0.2.0") < 0){
        stop(paste0("Please upgrade your autoCorrection installation to have ",
                "the required version. The current version is: ", aeVersion))
    }
    return(TRUE)
}

