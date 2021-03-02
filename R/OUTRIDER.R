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
OUTRIDER <- function(ods, q, controlData=TRUE, 
                    latent_space_model=c("autoencoder", "pca"), 
                    decoder_model=c("autoencoder", "pca"),
                    covariates=NULL, usePython=checkUsePython(ods, covariates), 
                    prepro_options=getDefaultPreproParams(ods),
                    pvalue_options=getDefaultPvalueParams(ods),
                    useBasilisk=FALSE,
                    implementation='autoencoder', 
                    BPPARAM=bpparam(), ...){
    checkOutrider2DataSet(ods)
    
    # pass on to the approriate implementation
    if(!missing(implementation)){
        warning("Using implementation argument is deprecated. Use 
                latent_space_model and decoder_model instead.")
        latent_space_model <- implementation
        decoder_model <- implementation
    }
    
    latent_space_model <- tolower(latent_space_model)
    decoder_model <- tolower(decoder_model)
    
    message(date(), ": Preprocessing ...")
    ods <- preprocess(ods, prepro_options=prepro_options)
    
    if(isTRUE(controlData)){
        message(date(), ": Controlling for confounders ...")
        ods <- controlForConfounders(ods, q=q, 
                                    prepro_options=prepro_options,
                                    latent_space_model=latent_space_model, 
                                    decoder_model=decoder_model,
                                    usePython=usePython, BPPARAM=BPPARAM, 
                                    covariates=covariates, 
                                    useBasilisk=useBasilisk, ...)
    }
    
    if(profile(ods) == "outrider" && 
            (isFALSE(controlData) | grepl("^(peer|pca)$", latent_space_model))){
        message(date(), ": Fitting NB to the data ...")
        ods <- fit(ods, BPPARAM=BPPARAM)
    }
    
    message(date(), ": P-value calculation ...")
    ods <- computePvalues(ods, BPPARAM=BPPARAM, 
                            distribution=pvalue_options$distribution,
                            method=pvalue_options$fdr_method, 
                            alternative=pvalue_options$pvalue_alternative)
    
    message(date(), ": Fold-change / zScore / delta calculation ...")
    ods <- computeEffectSizes(ods, distribution=pvalue_options$distribution,
                            effect_types=pvalue_options$effect_types, 
                            peerResiduals=grepl('^peer$', latent_space_model))
    
    validObject(ods)
    return(ods)
}

#' Preprocessing
#' 
#' Preprocesses an Outrider2DataSet according to the supplied parameter 
#' settings 
#' 
#' @param ods Outrider2DataSet
#' @param prepro_options List of preprocessing options. The default settings 
#'     can be obtained with getDefaultPreproParams(ods), and can be modified. 
#'     Recognized options for this function are 'prepro_func' (preprocessing 
#'     function, either NULL or a function), 'sf_norm' (TRUE or FALSE) and 
#'     'data_trans' (transformation function, one of 'log', 'log1p', or NULL).
#' @return An OutriderDataSet, with preprocessing output added.
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' prepro_opts <- getDefaultPreproParams(ods)
#' prepro_opts
#' ods <- preprocess(ods, prepro_options=prepro_opts)
#' 
#' @export
preprocess <- function(ods, prepro_options=getDefaultPreproParams(ods)){
    
    # apply preprocessing function if requested
    if(!is.null(prepro_options$prepro_func)){
        message("\t", date(), ":  Applying preprocessing function ...")
        prepro_func <- match.fun(prepro_options$prepro_func)
        preprocessed(ods) <- prepro_func(observed(ods))
    } else{
        preprocessed(ods) <- observed(ods)
    }
    
    # compute sizefactors and sf normalize if requested
    if(isTRUE(prepro_options$sf_norm)){
        message("\t", date(), ": SizeFactor estimation ...")
        ods <- estimateSizeFactors(ods)
    }
    else{
        sizeFactors(ods) <- rep(1, ncol(ods))
    }
    X_sf_norm <- t(t(preprocessed(ods)) / sizeFactors(ods))
    
    # compute transformed values
    if(!is.null(prepro_options$data_trans)){
        trans_func <- match.fun(prepro_options$data_trans)
        transformed(ods) <- trans_func(X_sf_norm)
    } else{
        transformed(ods) <- X_sf_norm
    }
    
    # add used options to metadata(ods)
    metadata(ods)$prepro_options <- prepro_options
    
    return(ods)
}

