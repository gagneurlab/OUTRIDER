autoCorrectED <- function(ods, q, theta=25, control=list(), debug=FALSE, ...){
    
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
    w_guess <- c(as.vector(pc),as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
               lossED(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    fit <- autoCorrectFit(w_guess, lossED, lossGradED, k, x, s, xbar, theta, 
                          control, ...)
    
    #Check that fit converged
    if(fit$convergence!=0){
        warning(paste0("Fit didn't converge with warning: ", fit$message))
    }
    
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               lossED(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictED(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}



lossED <- function(w, k, x, s, xbar, theta, ...){
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    ## encoding 
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    E <- W[,1:(ncol(W)/2)]
    D <- W[,(ncol(W)/2+1):ncol(W)]
    #theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    y <- t(t(x%*%E %*% t(D)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=theta, log=TRUE)
    - mean( ll )
}


lossGradED <- function(w, k, x, s, xbar, theta, ...){
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    E <- W[,1:(ncol(W)/2)]
    D <- W[,(ncol(W)/2+1):ncol(W)]
    theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    # dW:
    t1 <- t(x) %*% (k %*% D)
    #t1 <- armaMatMultAtBC(x, k, W)
    t2 <- t(k) %*% (x %*% E)
    #t2 <- armaMatMultAtBC(k, x, W)
    y <- t(t(x%*%E %*% t(D)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)
    kt <- (k + theta)*y_exp/(y_exp+theta)
    
    t3 <- t(x) %*% (kt %*% D)
    t4 <- t(kt) %*% (x %*% E)
    
    dE <- (-t1  + t3)/prod(dim(k))
    dD <- (- t2 + t4)/prod(dim(k))
    #db:
    db <- colSums(kt-k)/prod(dim(k))
    
    if(!all(is.finite(dE)) | !all(is.finite(dD))){
        browser()
    }
    return(c(dE,dD, db))
}


predictED <- function(w, k, s, xbar){
    x <-  t(t(log((1+k)/s)) - xbar)
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    E <- W[,1:(ncol(W)/2)]
    D <- W[,(ncol(W)/2+1):ncol(W)]
    
    y <- t(t(x%*%E %*% t(D)) +b + xbar)
    y_exp <- s*exp(y)
    y_exp
}