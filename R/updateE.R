
updateE <- function(ods, theta, control, BPPARAM, ...){
    e <- as.vector(getE(ods))
    D <- getD(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    x <- getx(ods)
    b <- getb(ods)
    
    fit <- optim(c(e,b), fn=lossE, gr=lossGradE, k=k, x=x, sf=sf, D=D,
                 theta=theta, method="L-BFGS-B", control=control, ...)
    
    # Check that fit converged
    if(fit$convergence!=0){
        warning(paste0("Fit didn't converge with warning: ", fit$message))
    }
    
    # update ods
    ods <- setE(ods, fit$par[1:length(e)])
    ods <- setb(ods, fit$par[1:length(b) + length(e)])
    
    return(ods)
}

lossE <- function(e, D, k, x, sf, theta, minMu=0.01, ...){
    
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    
    ## encoding 
    E <- getWeights(e, nr=ncol(k))
    b <- getBias(e, nr=ncol(k))
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- minMu + sf * exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=theta, log=TRUE)
    - mean( ll )
}


lossGradE <- function(e, D, k, x, sf, theta, minMu=0.01, ...){
    E <- getWeights(e, nr=ncol(k))
    b <- getBias(e, nr=ncol(k))
    theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    # dW:
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- minMu + sf * exp(y)
    kt <- (k + theta) * y_exp / (y_exp + theta)
    t1 <- t(x) %*% (k %*% D)
    t3 <- t(x) %*% (kt %*% D)
    
    # answers (dE and db)
    dE <- (-t1  + t3)/prod(dim(k))
    db <- colSums(kt-k)/prod(dim(k))
    
    return(c(dE, db))
}
