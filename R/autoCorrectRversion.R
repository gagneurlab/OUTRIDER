#' 
#' Autoencoder function to correct for confounders.
#' 
#' This is the wrapper function for the autoencoder implementation. 
#' It can be used to call the standard R implementation or the experimental
#' Python implementation.
#'
#' @param ods An OutriderDataSet object
#' @param q The encoding dimensions
#' @param theta The dispersion parameter
#' @param implementation "R", the default, will use the R implementation or 
#'             "python" to use the python/tensorflow experimental implementation
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} instance
#'             to be used for parallel computing.
#' @param ... passed on to the autoencoder implementing method. In the case of 
#'             the R implementation it is passed to the optim function. 
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' \dontshow{
#'     ods <- ods[1:10,1:10]
#' }
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
#' @export
autoCorrect <- function(ods, q, theta=25, 
                    implementation=names(autoEncoderImplList),
                    BPPARAM=bpparam(), ...){
    
    # error checking
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    checkSizeFactors(ods)
    
    if(!missing(q)){
        if(!is.numeric(q) && q > 1){
            stop("Please provide an integer greater then 1 for q.")
        }
        if(q >= nrow(ods)){
            stop("Please use a q smaller than the number of features.")
        }
        if(q >= ncol(ods)){
            stop("Please use a q smaller than the number of samples.")
        }
    } else {
        q <- getBestQ(ods)
        if(is.na(q)){
            q <- estimateBestQ(ods)
            message('Using the default for q with: ', q)
        }
    }
    
    # pass on to the correct implementation
    impl <- match.arg(implementation)
    if(is.null(impl) || is.na(impl)){
        stop("Requested autoCorrect implementation is unknown.")
    }
    message(date(), ": Using the ", impl, " implementation for autoCorrect.")
    aeFun <- autoEncoderImplList[[impl]]
    ans <- aeFun(ods=ods, q=q, theta=theta, BPPARAM=BPPARAM, ...)
    message(date(), ": Used the ", impl, " implementation for autoCorrect.")
    
    return(ans)
}

#' 
#' Autoencoder function to correct for confounders.
#'
#' @param ods An uormalized OUTRIDER data set
#' @param q the encoding dimension used.
#' @param theta value used in the likelihood (default=25).
#' 
#' @noRd
autoCorrectR <- function(ods, q, theta=25, control=list(), debug=FALSE, ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }

    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
                loss(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    fit <- autoCorrectFit(w_guess, loss, lossGrad, k, x, s, xbar, theta, 
            control, ...)
    
    #Check that fit converged
    if(fit$convergence!=0){
        warning(paste0("Fit didn't converge with warning: ", fit$message))
    }
    
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
                loss(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}

autoCorrectFit <- function(w, loss, lossGrad, k, x, s, xbar, theta, control, 
                           PCAonError=FALSE, ...){
    if(isTRUE(PCAonError)){
        fit <- NULL
        tryCatch({
            fit <- optim(w, loss, gr=lossGrad, k=k, x=x, s=s, xbar=xbar, 
                    theta=theta, method="L-BFGS-B", control=control, ...)
            },
            error = function(e) warning("Catched error: ", e$message))
        
        if(is.null(fit)){
            warning('An error occured during the autoencoder fit. ', 
                    'The initial PCA values are used.')
            fit <- list(
                convergence = 255,
                par = w,
                message = 'Errored during autoCorrect fitting.')
        }
        return(fit)
    }else{
        fit <- optim(w, loss, gr=lossGrad, k=k, x=x, s=s, xbar=xbar, 
                     theta=theta, method="L-BFGS-B", control=control, ...)
        return(fit)
    }
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
#' ods <- autoCorrect(ods)
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
