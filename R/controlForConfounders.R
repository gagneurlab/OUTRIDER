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
#' @export
controlForConfounders <- function(ods, q,
                    implementation=c('autoencoder', 'pca'), 
                    BPPARAM=bpparam(), ...){
    
    # error checking
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    checkSizeFactors(ods)
    
    if(!missing(q)){
        if(!(is.numeric(q) & q > 1 & q <= min(dim(ods)))){
            stop("Please provide for q an integer greater than 1 and smaller ", 
                    "than number of samples or genes.")
        }
    } else {
        q <- getBestQ(ods)
        if(is.na(q)){
            q <- estimateBestQ(ods)
            message('Using estimated q with: ', q)
        } else {
            message('Using provided q with: ', q)
        }
    }
    
    # pass on to the approriate implementation
    implementation <- tolower(implementation)
    implementation <- match.arg(implementation)
    aeFun <- switch(implementation,
        autoencoder = { 
            function(ods, q, ...){ fitAutoencoder(ods, q, ...) } },
        pca         = { 
            function(ods, q, BPPARAM, ...){ autoCorrectPCA(ods, q, ...) } },
        stop("Requested control implementation is unknown.")
    )
    
    message(date(), ": Using the ", implementation, 
            " implementation for controlling.")
    ans <- aeFun(ods=ods, q=q, BPPARAM=BPPARAM, ...)
    message(date(), ": Used the ", implementation, 
            " implementation for controlling.")
    
    return(ans)
}

#' OUTRIDER2 function for controlling for confounders, including option to 
#' fit in python.
#' 
#' @export
controlForConfounders2 <- function(ods, q,
                            implementation=c('autoencoder', 'pca'), 
                            usePython=FALSE, BPPARAM=bpparam(), ...){
    
    # error checking
    checkOutrider2DataSet(ods)
    checkDataRequirements(ods)
    checkPreprocessing(ods)
    checkFitInR(ods, !usePython)
    
    if(!missing(q)){
        if(!(is.numeric(q) & q > 1 & q <= min(dim(ods)))){
            stop("Please provide for q an integer greater than 1 and smaller ", 
                 "than number of samples or genes.")
        }
    } else {
        q <- getBestQ(ods)
        if(is.na(q)){
            q <- estimateBestQ(ods)
            message('Using estimated q with: ', q)
        } else {
            message('Using provided q with: ', q)
        }
    }
    
    # pass on to the approriate implementation
    implementation <- tolower(implementation)
    implementation <- match.arg(implementation)
    if(isTRUE(usePython)){
        aeFun <- function(ods, q, ...){ runPyCorrection(ods, q, 
                                            ncpus=bpnworkers(BPPARAM), ...) }
    } else{
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
    message(date(), ": Used the ", implementation, 
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
    stopifnot(is(ods, 'OutriderDataSet'))
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
runPyCorrection <- function(ods, q, covariates=NULL, 
                            ncpus=bpnworkers(bpparam()), ...){
    
    ### create args_param needed for py_outrider
    if(is(ods, "OutriderDataSet")){
        argsParam <- c("-p=outrider")
    }
    else if(is(ods, "ProtriderDataSet")){
        argsParam <- c("-p=protrider")
    }
    else{
        
        argsParam <- c(paste0("-dis=", modelParams(ods, "distribution")),
                        paste0("-dt=", modelParams(ods, "transformation")),
                        paste0("-pre=", modelParams(ods, "preprocessing")))    
    }
    
    argsParam <- c(argsParam, paste0("-q=", q), paste0("-cpu=", ncpus))
    
    if(!is.null(covariates)){
        argsParam <- c(argsParam,  paste0("-cov=", covariates))
    }
    
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
        # if("loss_distribution" %in% names(dots)){
        #     argsParam <- c(argsParam,  paste0("-ld=", dots$loss_distribution))
        # }
    }
    
    ### run python backend
    pyRes <- py_outrider$main$run_from_R_OUTRIDER(as.data.frame(t(raw(ods))), 
                                                  as.data.frame(colData(ods)), 
                                                  argsParam)
    # return(pyRes)
    
    ### add information from fit back to ods
    ods <- addPyFitResults(ods, pyRes)

    # return ods
    return(ods)
    
}

#' Annotates ods object with informations from python fit
#' @noRd
addPyFitResults <- function(ods, pyRes){
    
    # add fit information from python result back to ods
    normalizationFactors(ods) <- t(extract_pyRes(pyRes, "X_pred"))
    theta(ods) <- extract_pyRes(pyRes, "par_meas")
    
    # extract encoding dim (possibly with covariates)
    # q_with_cov <- extract_pyRes(pyRes, "encod_dim_cov")
    # 
    # # add encoder and decoder weights
    # E(ods) <- extract_pyRes(pyRes, "encoder_weights")
    # D(ods) <- extract_pyRes(pyRes, "decoder_weights")
    # b(ods) <- extract_pyRes(pyRes, "decoder_bias")
    
    # return ods 
    return(ods)
}

#' Function to extract a matrix or vector from the python result (xarray)
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

