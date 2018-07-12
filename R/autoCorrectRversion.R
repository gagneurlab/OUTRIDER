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
#' @param ... passed on to the autoencoder implementing method. In the case of 
#'             the R implementation it is passed to the optim function. 
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' \dontshow{
#'     ods <- ods[10:35,35:70]
#' }
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
#' @export
autoCorrect <- function(ods, q, theta=25, 
                    implementation=c("R", "python", "PEER"), ...){
    
    # error checking
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    checkSizeFactors(ods)
    
    if(!missing(q)){
        if(!is.numeric(q) && q > 0){
            stop("Please provide an integer greater then 0 for q.")
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
            q <- 5
        }
    }
    
    # pass on to the correct implementation
    if(match.arg(implementation)=='R'){
        return(autoCorrectR(ods, q, theta, ...))
    }
    if(match.arg(implementation)=='PEER'){
        return(peer(ods))
    }
    return(autoCorrectPython(ods, ...))
}

#' 
#' Autoencoder function to correct for confounders.
#'
#' @param ods An uormalized OUTRIDER data set
#' @param q the encoding dimension used.
#' @param theta value used in the likelihood (default=25).
#' 
#' @noRd
autoCorrectR <- function(ods, q, theta=25, control=list(), ...){
    
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
    fit <- optim(w_guess, loss, gr = lossGrad, k=k, x=x, s=s, xbar=xbar, 
            theta=theta, method="L-BFGS-B", control=control, ...)
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
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' computeLatentSpace(ods)[,1:6]
#' 
#' @export
computeLatentSpace <- function(ods){
    stopifnot(is(ods, 'OutriderDataSet'))
    if(!'weights' %in% names(metadata(ods))){
        stop('No weights are stored in this OutriderDataSet object. ',
                'Please compute weights before extracting the latent space.')
    }
    if(any(metadata(ods)[['dim']]!=dim(ods))){
        stop('The OutriderDataSet dimension changed and does not match with ',
                'the existing weights. Please recompute the weights. ',
                'Computation not possible on this data.')
    }
    
    # get data
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    w <- metadata(ods)[['weights']]
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    l <- t(x%*%W) 
    
    if(ncol(l)!=ncol(ods)){
        stop('Dimensions do not match.')
    }
    return(l)
}



## loss and gradient of loss function accesed by the functions above. ##

#' loss function
#'
#' @param w weight matrix 
#' @param k counts
#' @param x log-centered counts
#' @param s size factors
#' @param xbar offset 
#' @param theta value used in the likelihood (default=25).
#'
#' @return Returns the negative log likelihood of the negative binomial 
#' @noRd
loss <- function(w, k, x, s, xbar, theta){
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    
    b <- b + xbar
    truncLogLiklihood(k, x, W, b, s, theta)
}


#' gradient of loss function
#'
#' @param w weight matrix 
#' @param k counts
#' @param x log-centered counts
#' @param s size factors
#' @param xbar offset 
#' @param theta value used in the likelihood (default=25).
#'
#' @return returns the gradient of the loss function.
#' @noRd
lossGrad <- function(w, k, x, s, xbar, theta){
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    
    b <- b + xbar
    gradLogLiklihood(k, x, W, b, s, theta)
}



#'
#' Predict controlled Counts.
#' 
#' @param w weight matrix 
#' @param k counts
#' @param s size factors
#' @param xbar offset 
#'
#' @return Returns the predicted corrections (predicted means).
#' @noRd
predictC <- function(w, k, s, xbar){
    x <-  t(t(log((1+k)/s)) - xbar)
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    
    y <- t(t(x%*%W %*% t(W)) +b + xbar)
    s*exp(y)
}

#'
#' Get the bias term from the weights vector
#' 
#' @noRd
getBias <- function(w, nr){
    w[seq_len(nr)+(length(w)/nr-1)*nr]
}

#'
#' Get the weights matrix from the weight vector without the bias term
#' 
#' @noRd
getWeights <- function(w, nr){
    matrix(w, nrow=nr)[,seq_len(length(w)/nr-1)]
}

