autoCorrectRCooksIter3 <- function(ods, q, theta=25, control=list(), 
                                   BPPARAM=bpparam(), ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    
    k_no <-replaceOutliersCooksOutrider(k, q=1, BPPARAM=BPPARAM, ...)
    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
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
    
    w_fit <- w_guess
    for(i in 1:10){
        
        k_no <-replaceOutliersCooksOutrider(k,predictC(w_fit, k, s, xbar), q+1, 
                                            BPPARAM=BPPARAM, ...)
        x0 <- log((1+k_no)/s)
        x <- t(t(x0) - xbar)
        
        control$maxit <- 10    
        fit <- optim(w_fit, lossNonOutlier, gr = lossGradNonOutlier, k=k_no, x=x, s=s, xbar=xbar, 
                     theta=theta, method="L-BFGS-B", control=control, ...)
        
        w_fit <- fit$par
        message('Iteration ', i, ' loss: ', loss(w_fit, k, x, s, xbar, theta))
        
        #Check that fit converged
        if(fit$convergence!=0){
            warning(paste0("Fit didn't converge with warning: ", fit$message))
        }
        
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


lossNonOutlier <- function(w, k, x, s, xbar, theta, outlier){
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
    ll <- dnbinom(k, mu=y_exp, size=theta, log=TRUE)
    - sum( ll[!outlier] )
}




lossGradNonOutlier <- function(w, k, x, s, xbar, theta, outlier){
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    
    # dW:
    y <- t(t(x%*%W %*% t(W)) + xbar + b)
    y_exp <- s*exp(y)
    
    #k1 <- k*y_exp /pmax(y_exp-1,1E-10)
    k1 <- k
    k1[outlier] <- 0 
    t1 <- t(x) %*% (k1 %*% W)

    t2 <- t(k1) %*% (x %*% W)
    kt <- (k + theta)*(y_exp)/(y_exp + theta)
    #kt <- (k + theta)*(y_exp)/(pmax(y_exp-1,1E-10) + theta)
    kt[outlier] <- 0 
    t3 <- t(x) %*% (kt %*% W)
    #t3 <- armaMatMultAtBC(x, kt, W)
    t4 <- t(kt) %*% (x %*% W)
    #t4 <- armaMatMultAtBC(kt, x, W)
    dw <- (-t1 - t2 + t3 + t4)#/prod(dim(k))
    
    #db:
    db <- colSums(kt-k1)#/prod(dim(k))
    
    return(c(dw, db))
}
