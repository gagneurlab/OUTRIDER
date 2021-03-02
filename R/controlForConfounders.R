#' 
#' Autoencoder function to correct for confounders.
#' 
#' This is the wrapper function for the autoencoder implementation. 
#' It can be used to call the standard R or the python implementation.
#'
#' @param ods An OutriderDataSet object
#' @param q The encoding dimensions
#' @param prepro_options A list specifying the preprocessing options that 
#'          should be used. If the user wants to change some options, please 
#'          obtain the default list using \code{getDefaultPreproParams(ods)} 
#'          and adapt.
#' @param latent_space_model Specifies which latent space fitting model to use.
#'             "autoencoder", the default, will use the autoencoder
#'             implementation. Also 'pca' can be used to control
#'             for confounding effects
#' @param decoder_model Specifies which decoder fitting model to use.
#'             "autoencoder", the default, will use the autoencoder
#'             implementation. Also 'pca' can be used to control
#'             for confounding effects.
#' @param implementation Deprecated. Use latent_space_model and decoder_model 
#'             instead.
#' @param covariates Character vector indicating the known covariates that 
#'          should be considered explicitly in the fit. Defaults to NULL, 
#'          meaning no known covariates are used. Covariates given here must 
#'          be columns in \code{colData(ods)}. Only considered if 
#'          \code{usePython=TRUE}.
#' @param usePython Indicates whether the python or the R implementation of the
#'             given confounding method should be used. 
#' @param useBasilisk Specifies whether a conda environment installed with 
#'         basilik should be used to run the python code. If FALSE, 
#'         it is assumed that py_outrider is installed and the correct python 
#'         binary is specified using either 
#'         reticulate::use_python(..., force=TRUE) or 
#'         reticulate::use_condaenv(..., force=TRUE), so that py_outrider can 
#'         be loaded with reticulate::import("py_outrider").
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
                            prepro_options=getDefaultPreproParams(ods),
                            latent_space_model=c("autoencoder", "pca"), 
                            decoder_model=c("autoencoder", "pca"),
                            covariates=NULL, 
                            usePython=checkUsePython(ods, covariates),
                            useBasilisk=FALSE, 
                            implementation=latent_space_model, 
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
    if(!missing(implementation)){
        warning("Using implementation argument is deprecated. Use 
                latent_space_model and decoder_model instead.")
        latent_space_model <- implementation
        decoder_model <- latent_space_model
    }
    latent_space_model <- tolower(latent_space_model)
    latent_space_model <- match.arg(latent_space_model)
    metadata(ods)[["latent_space_model"]] <- latent_space_model
    decoder_model <- tolower(decoder_model)
    decoder_model <- match.arg(decoder_model)
    metadata(ods)[["decoder_model"]] <- decoder_model
    if(isTRUE(usePython)){
        aeFun <- function(ods, q, ...){ 
            runPyCorrection(ods, q, prepro_options,
                            latent_space_model=latent_space_model,
                            decoder_model=decoder_model, 
                            covariates=covariates,
                            useBasilisk=useBasilisk,
                            ncpus=bpnworkers(BPPARAM), ...) }
    } else{
        checkFitInR(ods, ...)
        if(latent_space_model != decoder_model){
            warning("The R implementation does not support different models\n",
            "for fitting the latent space and the decoder. Use the\n",
            "python version if you want do this. Running the R fit with\n",
            "latent_space_model = decoder_model = ", latent_space_model)
        }
        decoder_model <- latent_space_model
        aeFun <- switch(latent_space_model,
            autoencoder = { 
                function(ods, q, ...){ fitAutoencoder(ods, q, ...) } },
            pca         = { 
                function(ods, q, BPPARAM, ...){ autoCorrectPCA(ods, q, ...) } },
            stop("Requested control implementation is unknown.")
        )
    }
    
    implementation_string <- 
        paste(unique(c(latent_space_model, decoder_model)), collapse = " + ")
    progLang <- ifelse(isTRUE(usePython), "python", "R")
    message(date(), ": Using the ", progLang, " ", implementation_string, 
            " implementation for controlling.")
    ans <- aeFun(ods=ods, q=q, BPPARAM=BPPARAM, ...)
    message(date(), ": Used the ", progLang, " ", implementation_string, 
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
#' @param latent_space_model The model for fitting the latent space.
#' @param decoder_model The model for fitting the decoder.
#' @param prepro_options A list specifying the preprocessing options that 
#'          should be used. If the user wants to change some options, please 
#'          obtain the default list using \code{getDefaultPreproParams(ods)} 
#'          and adapt.
#' @param ncpus The number of cores for running the model fit in parallel.
#' @param ... Additional parameters passed onto the python implementation. 
#'          Currently supported are loss_distribution, iterations, convergence,
#'          seed and verbose.
#' 
#' @return An ods object including the control factors 
#' 
#' @noRd
runPyCorrection <- function(ods, q=NA, covariates=NULL, 
                            latent_space_model=c("autoencoder", "pca"), 
                            decoder_model=c("autoencoder", "pca"), 
                            prepro_options=getDefaultPreproParams(ods),
                            ncpus=bpnworkers(bpparam()), 
                            useBasilisk=FALSE, ...){
    
    # create parameter list needed for py_outrider run
    profile <- profile(ods)
    profile <- ifelse(profile == "other", "pca", profile)
    argsParam <- c(paste0("-p=", profile))
    
    # set correct input matrix
    prepro_func <- prepro_options$prepro_func
    if(!is.null(prepro_func)){
        inputMatrix <- as.data.frame(t(preprocessed(ods)))
        prepro_func <- "none"
    } else{
        inputMatrix <- as.data.frame(t(observed(ods)))
        prepro_func <- "none" 
    }
    
    # set latent_space fit model
    latent_space_model <- match.arg(latent_space_model)
    if(latent_space_model == "autoencoder"){
        argsParam <- c(argsParam, "--latent_space_model=AE")    
    } else if(latent_space_model == "pca"){
        argsParam <- c(argsParam, "--latent_space_model=PCA")    
    } else{
        stop("Unknown latent_space_model: ", latent_space_model)
    }
    
    # set decoder fit model
    decoder_model <- match.arg(decoder_model)
    if(decoder_model == "autoencoder"){
        argsParam <- c(argsParam, "--decoder_model=AE")    
    } else if(decoder_model == "pca"){
        argsParam <- c(argsParam, "--decoder_model=PCA")    
    } else{
        stop("Unknown decoder_model: ", decoder_model)
    }
    
    # set covariates to include in fit
    if(!is.null(covariates)){
        argsParam <- c(argsParam,  "-cov", covariates)
        covariates(ods) <- covariates
    }
    
    # check supplied data_trans
    data_trans <- prepro_options$data_trans
    if(!is.null(data_trans)){
        if(!(data_trans %in% c("log", "log1p"))){
            stop("cannot pass custom transform function to python. Please use 
                the python code directly or set data_trans to NULL, 'log' or 
                'log1p'.")
        }
    } else{
        data_trans <- "none"    
    }
    
    # set remaining options
    argsParam    <- c(argsParam, 
                        paste0("-pre=", prepro_func),
                        paste0("-dt=", data_trans),
                        paste0("-sf=", prepro_options$sf_norm),
                        paste0("-c=", "True"),
                        paste0("-nf=", prepro_options$noise_factor))
    
    dots <- list(...)
    if(length(dots) > 0){
        if("loss_distribution" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-ld=", dots$loss_distribution))
        }
        if("iterations" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-i=", dots$iterations))
        }
        if("verbose" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-v=", dots$verbose))
        }
        if("seed" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("-s=", dots$seed))
        }
        if("convergence" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("--convergence=", 
                                                dots$convergence))
        }
        if("batch_size" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("--batch_size=", 
                                                dots$batch_size))
        }
        if("nr_latent_space_features" %in% names(dots)){
            argsParam <- c(argsParam,  paste0("--nr_latent_space_features=", 
                                                dots$nr_latent_space_features))
        }
    }
    
    
    # set encoding dim. if q=NA, hyperParOpt in python code will be used
    if(!is.na(q)){
        argsParam <- c(argsParam, paste0("-q=", q), 
                        paste0("--num_cpus=", ncpus))
        extractHyperParOptResults <- FALSE
    } else{
        extractHyperParOptResults <- TRUE
        if(length(dots) > 0){
            if("convergence_hyper" %in% names(dots)){
                argsParam <- c(argsParam, paste0("--convergence_hyper=", 
                                                    dots$convergence_hyper))
            }
            if("max_iter_hyper" %in% names(dots)){
                argsParam <- c(argsParam, paste0("--max_iter_hyper=", 
                                                    dots$max_iter_hyper) )
            }
        }
    }
    
    ### run python backend
    if(isTRUE(useBasilisk)){
        message("Connecting to the py_outrider package using basilisk. ",
            "Set 'useBasilisk' = FALSE \nif you prefer to manually set the ",
            "python binary using 'reticulate'.")
        
        # proc <- basiliskStart(py_outrider_env)
        proc <- basiliskStart(outrider2_env)
        on.exit(basiliskStop(proc))
        ods <- basiliskRun(proc, function(ods, inputMatrix, argsParam, 
                                    extractHyperParOptResults){ 
            py_outrider <- import("py_outrider")
            pyRes <- py_outrider$main_run$run_from_R_OUTRIDER(
                inputMatrix,
                as.data.frame(colData(ods)),
                argsParam)
            return(addPyFitResults(ods, pyRes, extractHyperParOptResults))
        }, ods=ods, inputMatrix=inputMatrix, argsParam=argsParam, 
            extractHyperParOptResults=extractHyperParOptResults)
        
    } else{ # use current env loaded with reticulate

        message("Connecting to the py_outrider python package using reticulate",
            " (useBasilisk = FALSE)... \nIn case of errors, please make sure ",
            "to specify the right python binary when loading R.\n",
            "If you rather want us to automatically setup a conda ",
            "environment with 'py_outrider'\ninstalled using the 'basilisk'",
            " package, please use the argument 'useBasilisk = TRUE'.")
        
        py_outrider <- import("py_outrider")
        pyRes <- py_outrider$main_run$run_from_R_OUTRIDER(
            inputMatrix, 
            as.data.frame(colData(ods)), 
            argsParam)
        
        ### add information from fit back to ods
        ods <- addPyFitResults(ods, pyRes, extractHyperParOptResults)
    }

    # return ods
    return(ods)
    
}

#' Annotates ods object with informations from python fit
#' @noRd
addPyFitResults <- function(ods, pyRes, extractHyperParOptResults=FALSE){
    
    # add fit information from python result back to ods
    py_builtins <- import_builtins()
    layers <- py_builtins$dict(pyRes$layers)
    normalizationFactors(ods) <- t(layers[["X_predicted"]])
    
    # add theta for neg bin fits
    if("dispersions" %in% names(py_builtins$dict(pyRes$varm))){
        theta(ods) <- as.vector(pyRes$varm["dispersions"])
    }
    
    # add other info
    sizeFactors(ods) <- as.vector(pyRes$obsm["sizefactors"])
    preprocessed(ods, withDimnames=FALSE) <- t(layers[["X_prepro"]])
    if("X_transformed" %in% names(layers)){
        transformed(ods, withDimnames=FALSE) <- t(layers[["X_transformed"]])
    }
    
    if(isTRUE(extractHyperParOptResults)){
        encDimTable <- rbindlist(lapply(pyRes$uns[["hyperpar_table"]], 
                                        function(x){ as.data.table(x)}))
        setnames(encDimTable, "encod_dim", "encodingDimension")
        setnames(encDimTable, "prec_rec", "evaluationLoss")
        setnames(encDimTable, "loss", "fitLoss")
        encDimTable[,evalMethod:='aucPR']
        metadata(ods)[["encDimTable"]] <- encDimTable
        metadata(ods)[['optimalEncDim']] <- NULL
        metadata(ods)[['optimalEncDim']] <- getBestQ(ods)
    }
    
    # add encoder and decoder weights
    E(ods) <- pyRes$uns["E"]
    D(ods) <- pyRes$uns["D"]
    b(ods) <- pyRes$uns["bias"]
    
    # return ods 
    return(ods)
}

