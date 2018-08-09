

updateD <- function(ods, theta, control, BPPARAM, ...){
    d <- getD(ods)
    E <- getE(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    x <- getx(ods)
    b <- getb(ods)
    
    fit <- optim(c(d, b), fn=lossD, gr=lossGradD, k=k, x=x, sf=sf, E=E,
                 theta=theta, method="L-BFGS-B", control=control, ...)
    
    # Check that fit converged
    if(fit$convergence!=0){
        warning(paste0("Fit didn't converge with warning: ", fit$message))
    }
    
    # update ods
    ods <- setD(ods, fit$par[1:length(d)])
    ods <- setb(ods, fit$par[1:length(b) + length(d)])
    
    return(ods)
}



lossDOld <- function(d, E, k, x, sf, theta, minMu=0.01, ...){
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    
    ## encoding 
    D <- getWeights(d, nr=ncol(k))
    b <- getBias(d, nr=ncol(k))
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- minMu + sf * exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=theta, log=TRUE)
    - mean( ll )
}


lossGradDOld <- function(d, E, k, x, sf, theta, minMu=0.01, ...){
    D <- getWeights(d, nr=ncol(k))
    b <- getBias(d, nr=ncol(k))
    theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    # dW:
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- minMu + sf * exp(y)
    kt <- (k + theta) * y_exp / (y_exp + theta)
    t2 <- t(k) %*% (x %*% E)
    t4 <- t(kt) %*% (x %*% E)
    
    # answers (dD and db)
    dD <- (-t2 + t4)/prod(dim(k))
    db <- colSums(kt-k)/prod(dim(k))
    
    return(c(dD, db))
}


lossD <- function(d, k, H, sf, theta, minMu=0.01){
    b <- d[1]
    d <- d[-1]
    
    y <- H %*% d + b
    yexp <- minMu + sf * exp(y)
    
    ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    
    return(-ll)
}


gradD <- function(d, k, H, sf=1, theta, minMu=0.01){
    b <- d[1]
    d <- d[-1]
    
    y <- c(H %*% d + b)
    yexp <- minMu + sf * exp(y)
    t1 <- colMeans(c(k * sf * exp(y) / yexp) * H)
    
    kt <- (k + theta) * sf * exp(y) / (sf * exp(y) + theta)
    t2 <- colMeans(c(kt * sf * exp(y) / mu) * H) 
    
    dd <- t2-t1
    db <- mean(kt - k * sf * exp(y) / mu)
    
    return(c(db, dd))
}


debugLossD <- function(){
    samples <- 80
    q<- 2
    
    #Hidden Space is a sample time q matrix. 
    H <- matrix(c(rep(c(-1,1), each=samples/2), 
                  rep(c(-2,2), samples/2)), ncol=2)
    
    H
    D_true <- rnorm(q)
    y_true <- H %*% D_true + 3
    mu_true <- 0.01 + exp(y_true)
    
    k <- rnbinom(length(mu_true), mu = mu_true, size=25)
    
    #library(MASS)
    #glm.nb(k~H)
    
    init<-c(mean(log(k+1)),0,0)
    #b <- 3
    #d <- D_true
    
    
    fit <- optim(init, fn=lossD, gr=gradD, k=k, H=H, s=1, theta=25, method='L-BFGS')
    fit$par
    
    D_true
    lossD(k, init, H, 25)
    lossD(k,c(3, D_true), H,  25)
}

