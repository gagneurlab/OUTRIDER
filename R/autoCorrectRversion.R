sourceCpp("src/matMult.cpp")

#' 
#' Autoencoder function to correct for co-founders.
#' 
#' This is the wrapper function for the autoencoder implementation. 
#' It can be used to call the standard R implementation or the experimental
#' Python implementation.
#'
#' @param ods an OutriderDataSet object
#' @param q the encoding dimensions
#' @param theta the dispersion parameter
#' @param implementation "R", the default will use the R implementation or 
#'             "python" to use the python/tensorflow experimental implementation
#' @param ... passed on to the autoencoder implementing method.
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
#' @export
autoCorrect <- function(ods, q=20, theta=25, 
                    implementation=c("R", "python"), ...){
    
    # error checking
    if(!is(ods, 'OutriderDataSet')){
        stop('Please provide an OutriderDataSet')
    }
    if(q >= nrow(ods)){
        stop("Please use a q smaller than the number of features.")
    }
    if(q >= ncol(ods)){
        stop("Please use a q smaller than the number of samples.")
    }
    if(is.null(sizeFactors(ods))){
        stop(paste("Please calculate the size factors before calling", 
            "the autoCorrect function"))
    }
    
    # pass on to the correct implementation
    if(match.arg(implementation)=='R'){
        return(autoCorrectR(ods, q, theta))
    }
    return(autoCorrectPython(ods, ...))
}
    
#' 
#' Autoencoder function to correct for co-founders.
#'
#' @param ods An uormalized OUTRIDER data set
#' @param q the encoding dimension used.
#' @param theta value used in the likelihood (default=25).
#' 
#' @noRd
autoCorrectR <- function(ods, q=20, theta=25){

    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), rep(0,ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
            loss(w_guess, k, x, s, xbar, theta))
    )
    
    w_init <- c(as.vector(pc), rnorm(ncol(k), sd=0.001))
    # optimize log likelihood
    t <- Sys.time()
    fit <- optim(w_init, cmpLoss, gr = cmpLossGrad, k=k, x=x, s=s, xbar=xbar, 
        theta=theta, method="L-BFGS-B")
    w_fit <- fit$par
    print(paste0('Time elapsed: ', Sys.time() - t))
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


#' function to output the latentspace determined by the autoencoder.
#'
#' @param ods An OUTRIDER data set
#'
#' @return A matrix containing the by the autoencoder computed latent space.
#' @export
#'
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' computeLatentSpace(ods)
computeLatentSpace <- function(ods){
    stopifnot(is(ods, 'OutriderDataSet'))
    if(metadata(ods)[['dim']]!=dim(ods)){
        stop('The ods dimension changed. Computation not possible.')
    }
    # get data
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    weights <- metadata(ods)[['weights']]
    
    W <- matrix(weights, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
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
    ## log, size factored, and centered counts
    #x <-  t(t(log((1+k)/s)) - xbar)
    ## encoding
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]

    y <- t(t(x%*%W %*% t(W)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)

    ## log likelihood - truncated
    - mean(k * (log(s)+y)) + mean((k + theta) * log(y_exp + theta))
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
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    
    # dW:
    #t1 <- t(x) %*% k %*% tWe
    t1 <- armaMatMultAtBC(x, k, W)
    #t2 <- t(k) %*% x %*% tWe
    t2 <- armaMatMultAtBC(k, x, W)
    y <- t(t(x%*%W %*% t(W)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)
    kt <- (k + theta)*y_exp/(y_exp+theta)
    #t3 <- t(x) %*% kt %*% tWe
    t3 <- armaMatMultAtBC(x, kt, W)
    #t4 <- t(kt) %*% x %*% tWe
    t4 <- armaMatMultAtBC(kt, x, W)
    dw <- (-t1 - t2 + t3 + t4)/prod(dim(k))
    
    #db:
    db <- colSums(kt-k)/prod(dim(k))
    
    return(c(dw, db))
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
#' 
predictC <- function(w, k, s, xbar){
    x <-  t(t(log((1+k)/s)) - xbar)
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    y <- t(t(x%*%W %*% t(W)) +b + xbar)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    s*exp(y)
}


numericLossGrad <- function(fn, epsilon, w,...){
    grad <- list()
    for(i in seq_along(w)){
        eps <- rep(0, length(w))
        eps[i] <- epsilon
        grad[i] <- (fn(w + eps, ...) - fn(w -eps, ...))/2*epsilon
    }
    return(grad)
}

# Use the R compiler for the loss and grad functions.
cmpLossGrad <- cmpfun(lossGrad)
cmpLoss <- cmpfun(loss)

