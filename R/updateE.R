
updateE <- function(ods, theta, control, BPPARAM, ...){
    e <- as.vector(getE(ods))
    D <- getD(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    x <- getx(ods)
    b <- getb(ods)
    control$trace<-3
    
    fit <- optim(e, fn=truncLogLiklihoodE, gr=gradientE, k=k, x=x, sf=sf, D=D, 
            b=b, theta=theta, method="L-BFGS-B", control=control, ...)
    
    # Check that fit converged
    if(fit$convergence!=0){
        print(paste('updateE did not converge: ', fit$message))
        warning(paste0("Fit didn't converge with warning: ", fit$message))
    }
    
    # update ods
    ods <- setE(ods, fit$par)
    
    return(ods)
}

lossE <- function(e, D, k, b, x, sf, theta, minMu=0.01, ...){
    
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    
    ## encoding 
    E <-matrix(e, nrow=ncol(k))
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * (minMu + exp(y))
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=theta, log=TRUE)
    - mean( ll )
}

lossEtrunc <- function(e, D, k, b, x, sf, theta, minMu=0.01, ...){
    E <-matrix(e, nrow=ncol(k))
    
    y <- t(t(x %*% E %*% t(D)) + b)
    
    #ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    #ll = mean(k * log(yexp) - (k + theta)*log(yexp + theta))
    
    t1 <- k * (log(sf) + y + log(1 + minMu/exp(y)))
    t2 <- (k + theta) * (log(sf) + y + log(1 + minMu/exp(y))  + log(1+theta/(sf * (minMu + exp(y))))  )
    ll <- mean(t1 - t2)
    
    # if(!is.finite(ll) & debugMyCode){
    #   browser()
    # }
    
    return(-ll)
}


lossGradE <- function(e, D, k, b, x, sf, theta, minMu=0.01, ...){
    E <-matrix(e, nrow=ncol(k))
    theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    # dW:
    y <- t(t(x %*% E %*% t(D)) + b)

    yexp <- sf * (minMu + exp(y))
    #k1 <- k *sf* exp(y) / yexp         
    k1 <- k / (1 + minMu/exp(y) )          
    #kt <- (k + theta) * sf * exp(y) / (yexp + theta)
    kt <- (k + theta) / ( 1 + (minMu + theta/sf)/exp(y) )
    t1 <- t(x) %*% (k1 %*% D)
    t3 <- t(x) %*% (kt %*% D)
    
    # answers dE 
    dE <- (-t1 + t3)/prod(dim(k))
    
    return(dE)
}




##### Debug code
if(FALSE){
    numericLossGrad <- function(fn, epsilon, w,...){
        grad <- numeric(length(w))
        for(i in seq_along(w)){
            eps <- integer(length(w))
            eps[i] <- epsilon
            grad[i] <- (fn(w + eps, ...) - fn(w -eps, ...))/(2*epsilon)
        }
        return(grad)
    }
    plot(numericLossGrad(lossE, 1E-8, w=e, D=D, k=k, b=b, x=x, sf=sf, theta=theta, minMu=0.01),
         lossGradE(e, D, k, b, x, sf, theta, minMu=0.01));abline(0,1)
}




