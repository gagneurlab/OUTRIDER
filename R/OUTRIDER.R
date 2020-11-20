#' 
#' OUTRIDER - Finding expression outlier events
#' 
#' @description The OUTRIDER function runs the default OUTRIDER pipeline 
#' combinig the fit, the computation of Z scores and P-values.
#' All computed values are returned as an OutriderDataSet object.
#' 
#' To have more control over each analysis step, one can call each 
#' function separately.
#' 
#' \enumerate{
#'     \item \code{\link{estimateSizeFactors}} to calculate the sizeFactors
#'     \item \code{\link{controlForConfounders}} to control for 
#'               confounding effects
#'     \item \code{\link{fit}} to fit the negative binomial model 
#'               (only needed if the autoencoder is not used)
#'     \item \code{\link{computePvalues}} to calculate the nominal and 
#'               adjusted P-values
#'     \item \code{\link{computeZscores}} to calculate the Z scores
#' }
#' 
#' @inheritParams controlForConfounders
#' @param controlData If TRUE, the default, the raw counts are controled 
#'             for confounders by the autoencoder
#' @param ... Further arguments passed on to \code{controlForConfounders}
#' @return OutriderDataSet with all the computed values. The values are stored
#'             as assays and can be accessed by: \code{assay(ods, 'value')}.
#'             To get a full list of calculated values run:
#'             \code{assayNames(ods)}
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' implementation <- 'autoencoder'
#' \dontshow{
#'     ods <- ods[1:10,1:10]
#'     implementation <- 'pca'
#' }
#' ods <- OUTRIDER(ods, implementation=implementation)
#' 
#' pValue(ods)[1:10,1:10]
#' res <- results(ods, all=TRUE)
#' res
#' 
#' plotAberrantPerSample(ods)
#' plotVolcano(ods, 1)
#' 
#' @export
# OUTRIDER <- function(ods, q, controlData=TRUE, implementation='autoencoder', 
#                     BPPARAM=bpparam(), ...){
#     checkOutriderDataSet(ods)
#     implementation <- tolower(implementation)
#     
#     message(date(), ": SizeFactor estimation ...")
#     ods <- estimateSizeFactors(ods)
#     
#     if(isTRUE(controlData)){
#         message(date(), ": Controlling for confounders ...")
#         ods <- controlForConfounders(ods, q=q, 
#                 implementation=implementation, BPPARAM=BPPARAM, ...)
#     }
#     
#     if(isFALSE(controlData) | grepl("^(peer|pca)$", implementation)){
#         message(date(), ": Fitting the data ...")
#         ods <- fit(ods, BPPARAM=BPPARAM)
#     }
#     
#     message(date(), ": P-value calculation ...")
#     ods <- computePvalues(ods, BPPARAM=BPPARAM)
#     
#     message(date(), ": Zscore calculation ...")
#     ods <- computeZscores(ods, 
#             peerResiduals=grepl('^peer$', implementation))
#     
#     validObject(ods)
#     return(ods)
# }
# 
OUTRIDER <- function(ods, q, controlData=TRUE, implementation='autoencoder', 
                        usePython=ifelse(ncol(ods) < 500, FALSE, TRUE), 
                        BPPARAM=bpparam(), ...){
    checkOutrider2DataSet(ods)
    implementation <- tolower(implementation)
    
    message(date(), ": Preprocessing ...")
    ods <- preprocess(ods)
    
    if(isTRUE(controlData)){
        message(date(), ": Controlling for confounders ...")
        ods <- controlForConfounders(ods, q=q, implementation=implementation, 
                                        usePython=usePython, BPPARAM=BPPARAM, 
                                        ...)
    }
    
    if(modelParams(ods)$distribution == "negative binomial" && 
       (isFALSE(controlData) | grepl("^(peer|pca)$", implementation))){
        message(date(), ": Fitting NB to the data ...")
        ods <- as(ods, "OutriderDataSet")
        ods <- fit(ods, BPPARAM=BPPARAM)
    }
    
    message(date(), ": P-value calculation ...")
    ods <- computePvalues(ods, BPPARAM=BPPARAM)
    
    message(date(), ": Zscore and fold-change calculation ...")
    ods <- computeZscores(ods,
                          peerResiduals=grepl('^peer$', implementation))
    
    validObject(ods)
    return(ods)
}

#' @noRd
preprocess <- function(ods, normalized=FALSE){
    prepro <- modelParams(ods, "preprocessing")
    trans <- modelParams(ods, "transformation")
    sf_norm <- modelParams(ods, "sf_norm")
        
    # calculate sizefactors if needed
    if(isTRUE(sf_norm) || prepro == "vst"){
        message("\t", date(), ": SizeFactor estimation ...")
        ods <- estimateSizeFactors(ods)
    }
    
    # apply specified preprocessing
    if(prepro == "log"){
        message("\t", date(), ":  Taking the log of the data ...")
        if(sf_norm){
            sf <- sizeFactors(ods)
        } else{
            sf <- rep(1, ncol(ods))
        }
        preprocessed(ods) <- 
            log((observed(ods, normalized=normalized) + 1) /  sf)
    }
    if(prepro == "vst"){
        message("\t", date(), ":  Variance stabilizing using DESeq2 ...")
        preprocessed(ods) <- 
            DESeq2::varianceStabilizingTransformation(
                observed(ods, normalized=normalized))
    }
    
    return(ods)
}

