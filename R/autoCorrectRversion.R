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
#' @param ods 
#' @param q 
#' @param theta 
#' 
#' @noRd
autoCorrectR <- function(ods, q=20, theta=25){
    
    # get data
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) ## could also use ppca
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
    #pars <- cmpOptGD(w_guess, k, x, s, xbar, theta)
    fit <- optim(w_init, cmpLoss, gr = cmpLossGrad, k=k, x=x, s=s, xbar=xbar, 
                 theta=theta, method="L-BFGS-B")
    pars <- fit$par
    print(paste0('Time elapsed: ', Sys.time() - t))
    print(
        paste0('nb-PCA loss: ',
               loss(pars,k, x, s,xbar, theta))
    )

    correctionFactors <- t(predictC(pars, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- pars
    validObject(ods)
    return(ods)
}

computeLatentSpace <- function(ods){
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
        stop('Error')
    }
    
    return(l)
}


autoCorrectPCA <- function(ods, q=20, theta=25){
    # error checking
    if(is.null(sizeFactors(ods))){
        stop(paste("Please calculate the size factors before calling", 
                   "the autoCorrect function"))
    }
    # get data
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize using PCA
    pca <- pca(x, nPcs = q) ## could also use ppca
    pc  <- loadings(pca)
    w_guess <- as.vector(pc)
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
               loss(c(w_guess, rep(0,ncol(k))), k, x, s, xbar, theta))
    )
    
    pars <- c(w_guess, rep(0,ncol(k)))
    
    corrected <- predictC(pars, k, s, xbar)
    
    correctionFactors <- t(corrected)
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    validObject(ods)
    return(ods)
}


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

    ## log likelihood 
    ll <- dnbinom(k, mu= y_exp, size=theta, log=TRUE)
    - mean( ll )
}

predictC <- function(w, k, s, xbar){
    x <-  t(t(log((1+k)/s)) - xbar)
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    y <- t(t(x%*%W %*% t(W)) +b + xbar)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    s*exp(y)
}

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

numericLossGrad <- function(fn, epsilon, w,...){
    grad <- list()
    for(i in seq_along(w)){
        eps <- rep(0, length(w))
        eps[i] <- epsilon
        grad[i] <- (fn(w + eps, ...) - fn(w -eps, ...))/2*epsilon
    }
    return(grad)
}

# try to use the R compiler for the loss and grad functions.
cmpLossGrad <- cmpfun(lossGrad)
cmpLoss <- cmpfun(loss)

