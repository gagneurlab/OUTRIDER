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
#' @param subsets A named list of named lists specifying any number of gene 
#'             subsets (can differ per sample). For each subset, FDR correction
#'             will be limited to genes in the subset, and the FDR corrected 
#'             pvalues stored as an assay in the ods object in addition to the 
#'             transcriptome-wide FDR corrected pvalues. See the examples for 
#'             how to use this argument. 
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
#' # example of restricting FDR correction to subsets of genes of interest 
#' genesOfInterest <- list("sample_1"=sample(rownames(ods), 3), 
#'                          "sample_2"=sample(rownames(ods), 8), 
#'                          "sample_6"=sample(rownames(ods), 5))
#' genesOfInterest
#' ods <- OUTRIDER(ods, subsets=list("exampleSubset"=genesOfInterest))
#' padj(ods, subsetName="exampleSubset")[1:10,1:10]
#' res <- results(ods, all=TRUE)
#' res
#' 
#' @export
OUTRIDER <- function(ods, q, controlData=TRUE, implementation='autoencoder', 
                     subsets=NULL, BPPARAM=bpparam(), ...){
    checkOutriderDataSet(ods)
    implementation <- tolower(implementation)
    
    message(date(), ": SizeFactor estimation ...")
    ods <- estimateSizeFactors(ods)
    
    if(isTRUE(controlData)){
        message(date(), ": Controlling for confounders ...")
        ods <- controlForConfounders(ods, q=q, 
                implementation=implementation, BPPARAM=BPPARAM, ...)
    }
    
    if(isFALSE(controlData) | grepl("^(peer|pca)$", implementation)){
        message(date(), ": Fitting the data ...")
        ods <- fit(ods, BPPARAM=BPPARAM)
    }
    
    message(date(), ": P-value calculation ...")
    ods <- computePvalues(ods, subsets=subsets, BPPARAM=BPPARAM)
    
    message(date(), ": Zscore calculation ...")
    ods <- computeZscores(ods, 
            peerResiduals=grepl('^peer$', implementation))
    
    validObject(ods)
    return(ods)
}

