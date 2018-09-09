#' 
#' OUTRIDER - Finding aberrant expression events
#' 
#' @description The OUTRIDER function runs the default OUTRIDER pipline. 
#' Combinig the fit, the computation of z scores and P-values.
#' All computed values are returned as an OutriderDataSet object.
#' 
#' To have more control over each analysis step one can call each 
#' function seperatly.
#' 
#' \enumerate{
#'     \item \code{\link{estimateSizeFactors}} to calculte the sizeFactors
#'     \item \code{\link{autoCorrect}} to correct for confounding effects
#'     \item \code{\link{fit}} to fit the negative binomial model
#'     \item \code{\link{computePvalues}} to calculate the nominal and 
#'               adjusted P-values
#'     \item \code{\link{computeZscores}} to calculate the z scores
#' }
#' 
#' @param object An OutriderDataSet object containing the counts
#' @inheritParams autoCorrect
#' @param autoCorrect If TRUE, the default, the raw read counts are controled 
#'             for confounders by the autoencoder
#' @return OutriderDataSet with all the computed values. The values are stored
#'             as assays and can be accessed by: \code{assays(ods)[['value']]}.
#'             To get a full list of calculated values run:
#'             \code{assayNames(ods)}
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' runAutoCorrect <- TRUE
#' \dontshow{
#'     ods <- ods[1:10,1:10]
#'     runAutoCorrect <- FALSE
#' }
#' ods <- OUTRIDER(ods, autoCorrect=runAutoCorrect)
#' 
#' assays(ods)[['pValue']][1:10,1:10]
#' res <- results(ods, all=TRUE)
#' res
#' 
#' plotAberrantPerSample(ods)
#' plotVolcano(ods, 1)
#' 
#' @export
OUTRIDER <- function(object, q, autoCorrect=TRUE, implementation='edNew', 
                    BPPARAM=bpparam(), ...){
    implementation <- tolower(implementation)
    
    message(date(), ": SizeFactor estimation ...")
    object <- estimateSizeFactors(object)
    
    if(isTRUE(autoCorrect)){
        message(date(), ": Running autoencoder ...")
        object <- autoCorrect(object, q=q, implementation=implementation, 
                BPPARAM=BPPARAM, ...)
    }
    
    if(isFALSE(autoCorrect) | grepl("^(peer|pca)$", implementation)){
        message(date(), ": Fitting the data ...")
        object <- fit(object, BPPARAM=BPPARAM)
    }
    
    message(date(), ": P-value calculation ...")
    object <- computePvalues(object, BPPARAM=BPPARAM)
    
    message(date(), ": Zscore calculation ...")
    object <- computeZscores(object, 
            peerResidual=grepl('^peer$', implementation))
    
    validObject(object)
    return(object)
}

