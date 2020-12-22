#' 
#' Autoencoder function to correct for confounders.
#' 
#' This is the wrapper function for the autoencoder implementation. 
#' It can be used to call the standard R implementation or the experimental
#' Python implementation.
#'
#' @param ods An OutriderDataSet object
#' @param q The encoding dimensions
#' @param implementation "autoencoder", the default, will use the autoencoder
#'             implementation. Also 'pca' and 'peer' can be used to control
#'             for confounding effects
#' @param usePython Indicates whether the python or the R implementation of the
#'             given confounding method should be used. 
#' @param BPPARAM A 
#'     \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'             instance to be used for parallel computing.
#' @param ... Further arguments passed on to the specific implementation method.
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' implementation <- 'autoencoder'
#' \dontshow{
#'     ods <- ods[1:10,1:10]
#'     implementation <- 'pca'
#' }
#' ods <- estimateSizeFactors(ods)
#' ods <- controlForConfounders(ods, implementation=implementation)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
# #' @export
# controlForConfounders <- function(ods, q,
#                     implementation=c('autoencoder', 'pca'), 
#                     BPPARAM=bpparam(), ...){
#     
#     # error checking
#     checkOutriderDataSet(ods)
#     checkCountRequirements(ods)
#     checkSizeFactors(ods)
#     
#     if(!missing(q)){
#         if(!(is.numeric(q) & q > 1 & q <= min(dim(ods)))){
#             stop("Please provide for q an integer greater than 1 and smaller ", 
#                     "than number of samples or genes.")
#         }
#     } else {
#         q <- getBestQ(ods)
#         if(is.na(q)){
#             q <- estimateBestQ(ods)
#             message('Using estimated q with: ', q)
#         } else {
#             message('Using provided q with: ', q)
#         }
#     }
#     
#     # pass on to the approriate implementation
#     implementation <- tolower(implementation)
#     implementation <- match.arg(implementation)
#     aeFun <- switch(implementation,
#         autoencoder = { 
#             function(ods, q, ...){ fitAutoencoder(ods, q, ...) } },
#         pca         = { 
#             function(ods, q, BPPARAM, ...){ autoCorrectPCA(ods, q, ...) } },
#         stop("Requested control implementation is unknown.")
#     )
#     
#     message(date(), ": Using the ", implementation, 
#             " implementation for controlling.")
#     ans <- aeFun(ods=ods, q=q, BPPARAM=BPPARAM, ...)
#     message(date(), ": Used the ", implementation, 
#             " implementation for controlling.")
#     
#     return(ans)
# }

#' OUTRIDER2 function for controlling for confounders, including option to 
#' fit in python.
#' 
#' @export
controlForConfounders <- function(ods, q,
                            implementation=c('autoencoder', 'pca'), 
                            covariates=NULL, 
                            usePython=checkUsePython(ods, covariates),
                            BPPARAM=bpparam(), ...){
    
    # error checking
    checkOutrider2DataSet(ods)
    checkDataRequirements(ods)
    checkPreprocessing(ods)
    
    # check covariates argument
    if(!is.null(covariates) && isFALSE(usePython)){
        warning("Known covariates can only be included in the python ",
                "implementation. Setting usePython=TRUE.\n",
                "If you would like to use the R implementation, set ",
                "'covariates=NULL' in the function call.")
        usePython=TRUE
    }
    
    if(!missing(q) && !is.null(q)){
        if(!(is.numeric(q) & q > 1 & q <= min(dim(ods)))){
            stop("Please provide for q an integer greater than 1 and smaller ", 
                 "than number of samples or genes.")
        }
        message('Using specified q = ', q)
    } else {
        q <- getBestQ(ods)
        if(is.na(q)){
            if(isTRUE(usePython)){
                message("Estimating q with the python implementation")
            } else{
                q <- estimateBestQ(ods)
                message('Using estimated q with: ', q)
            }
        } else {
            message('Using provided q with: ', q)
        }
    }
    
    # pass on to the approriate implementation
    implementation <- tolower(implementation)
    implementation <- match.arg(implementation)
    metadata(ods)[["fitModel"]] <- implementation
    if(isTRUE(usePython)){
        aeFun <- function(ods, q, ...){ runPyCorrection(ods, q, implementation,
                                            covariates=covariates,
                                            ncpus=bpnworkers(BPPARAM), ...) }
    } else{
        checkFitInR(ods)
        aeFun <- switch(implementation,
            autoencoder = { 
                function(ods, q, ...){ fitAutoencoder(ods, q, ...) } },
            pca         = { 
                function(ods, q, BPPARAM, ...){ autoCorrectPCA(ods, q, ...) } },
            stop("Requested control implementation is unknown.")
        )
    }
    
    progLang <- ifelse(isTRUE(usePython), "python", "R")
    message(date(), ": Using the ", progLang, " ", implementation, 
            " implementation for controlling.")
    ans <- aeFun(ods=ods, q=q, BPPARAM=BPPARAM, ...)
    message(date(), ": Used the ", progLang, " ", implementation, 
            " implementation for controlling.")
    
    return(ans)
}

#' 
#' Extracting the latent space
#' 
#' Extracts the latent space from the OutriderDataSet object 
#' determined by the autoencoder.
#'
#' @param ods An OutriderDataSet
#'
#' @return A matrix containing the latent space determined by the autoencoder.
#'
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' \dontshow{
#'     ods <- ods[1:10, 1:10]
#' }
#' ods <- estimateSizeFactors(ods)
#' ods <- controlForConfounders(ods, implementation="pca")
#' computeLatentSpace(ods)[,1:6]
#' 
#' @export
computeLatentSpace <- function(ods){
    stopifnot(is(ods, 'Outrider2DataSet'))
    if(is.null(D(ods)) | is.null(E(ods))){
        stop('The D or E weights are not computed yet. Please fit the ',
                'autoencoder before extracting the latent space.')
    }
    if(any(metadata(ods)[['dim']] != dim(ods))){
        stop('The OutriderDataSet dimension changed and does not match with ',
                'the existing autoencoder fit. Please refit the autoencoder.')
    }
    
    return(H(ods))
}


#' 
#' Controlling for confounders using python backend
#' 
#' @param ods An OutriderDataSet
#' @param q The encoding dimension. 
#' @param covariates Character vector indicating the known covariates that 
#'          should be considered explicitly in the fit. Defaults to NULL, 
#'          meaning no known covariates are used. Covariates given here must 
#'          be columns in \code{colData(ods)}.
#' @param ncpus
#' 
#' @return An ods object including the control factors 
#' 
#' @noRd
runPyCorrection <- function(ods, q=NA, implementation=c("autoencoder", "pca"), 
                            covariates=NULL, ncpus=bpnworkers(bpparam()), ...){
    
    ### create parameter list needed for py_outrider run
    if(implementation == "autoencoder"){
        argsParam <- c("-p=outrider")    
    } else if(implementation == "pca"){
        argsParam <- c("-p=pca")    
    } else{
        stop("Unknown latent space fitting implementation: ", implementation)
    }

    distr <- switch(modelParams(ods, "distribution"),
                    "negative binomial"="neg_bin",
                    "gaussian"="gaus")
    sf_norm <- modelParams(ods, "sf_norm")        
    trans <- switch(modelParams(ods, "transformation"),
                    "log"  = ifelse(sf_norm, "sf", "log"),
                    "none" = "none")
    prepro <- switch(modelParams(ods, "preprocessing"),
                     "none" = "none",
                     "log"  = "sf_log",
                     "vst"  = "none")
    argsParam <- c(argsParam, paste0("-dis=", distr), paste0("-ld=", distr),
                    paste0("-dt=", trans),
                    paste0("-pre=", prepro))    
    
    # set encoding dim. if q=NA, hyperParOpt in python code will be used
    if(!is.na(q)){
        argsParam <- c(argsParam, paste0("-q=", q), paste0("--num_cpus=", ncpus))
        extractHyperParOptResults <- FALSE
    } else{
        extractHyperParOptResults <- TRUE
    }
    
    # set covariates to include in fit
    if(!is.null(covariates)){
        argsParam <- c(argsParam,  "-cov", covariates)
        covariates(ods) <- covariates
    }
    
    # set optional additional parameters
    dots <- list(...)
    if(length(dots) > 0){
        if("iterations" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-m=", dots$iterations))
        }
        if("verbose" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-v=", dots$verbose))
        }
        if("seed" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-s=", dots$seed))
        }
        if("loss_distribution" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-ld=", dots$loss_distribution))
        }
    }
    
    ### set correct input matrix
    if(modelParams(ods, "preprocessing") == "vst"){
        inputMatrix <- as.data.frame(t(preprocessed(ods)))
    } else{
        inputMatrix <- as.data.frame(t(observed(ods, normalized=FALSE)))
    }
    
    ### run python backend
    pyRes <- py_outrider$main$run_from_R_OUTRIDER(
        inputMatrix, 
        as.data.frame(colData(ods)), 
        argsParam)
    
    ### add information from fit back to ods
    ods <- addPyFitResults(ods, pyRes, extractHyperParOptResults)

    # return ods
    return(ods)
    
}

#' Annotates ods object with informations from python fit
#' @noRd
addPyFitResults <- function(ods, pyRes, extractHyperParOptResults=FALSE){
    
    # add fit information from python result back to ods
    normalizationFactors(ods) <- t(extract_pyRes(pyRes, "X_pred"))
    
    # add theta for neg bin fits
    if(modelParams(ods, "distribution") == "negative binomial"){
        theta(ods) <- extract_pyRes(pyRes, "par_meas")
    }
    
    if(isTRUE(extractHyperParOptResults)){
        encDimTable <- rbindlist(lapply(pyRes[["hyperpar_table"]], 
                                        function(x){ as.data.table(x)}))
        setnames(encDimTable, "encod_dim", "encodingDimension")
        setnames(encDimTable, "prec_rec", "evaluationLoss")
        setnames(encDimTable, "loss", "fitLoss")
        encDimTable[,evalMethod:='aucPR']
        metadata(ods)[["encDimTable"]] <- encDimTable
        metadata(ods)[['optimalEncDim']] <- NULL
        metadata(ods)[['optimalEncDim']] <- getBestQ(ods)
    }
    
    # extract encoding dim (possibly with covariates)
    # q_with_cov <- extract_pyRes(pyRes, "encod_dim_cov")
     
    # add encoder and decoder weights
    E(ods) <- extract_pyRes(pyRes, "encoder_weights")
    D(ods) <- extract_pyRes(pyRes, "decoder_weights")
    b(ods) <- extract_pyRes(pyRes, "decoder_bias")
    
    # return ods 
    return(ods)
}

#' Function to extract a matrix or vector from the py_outrider result (xarray)
#' @noRd
extract_pyRes <- function(pyRes, aname){
    assay_xr <- pyRes[[aname]]
    assay <- assay_xr$values
    
    if(length(dim(assay)) == 2){
        rownames(assay) <- assay_xr$coords$get("sample")$values
        colnames(assay) <- assay_xr$coords$get("meas")$values
    } else{
        names(assay) <- assay_xr$coords$get("meas")$values
    }
    return(assay)
}

