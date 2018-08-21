#'
#' Update E step for the autoencoder fit
#' 
#' @noRd
updateE <- function(ods, control, BPPARAM){
    e <- as.vector(E(ods))
    D <- D(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    x <- x(ods)
    b <- b(ods)
    theta <- theta(ods)
    mask <- t(exclusionMask(ods))
    
    control$trace <- 3
    fit <- optim(e, fn=truncLogLiklihoodE, gr=gradientE,
            k=k, x=x, sf=sf, D=D, b=b, theta=theta, exclusionMask=mask, 
            method="L-BFGS-B", lower=-100, upper=100, control=control)
    
    # Check that fit converged
    if(fit$convergence!=0){
        print(paste('Update E did not converge: ', fit$message))
    }
    
    # update ods
    E(ods) <- fit$par
    
    return(ods)
}

lossE <- function(e, D, k, b, x, sf, theta, ...){
    
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    
    ## encoding 
    E <-matrix(e, nrow=ncol(k))
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=theta, log=TRUE)
    - mean( ll )
}

lossEtrunc <- function(e, D, k, b, x, sf, theta, minMu=0, ...){
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


lossEtruncNonOutlier <- function(e, D, k, b, x, sf, theta, minMu=0, exlusionMask, ...){
    E <-matrix(e, nrow=ncol(k))
  
    y <- t(t(x %*% E %*% t(D)) + b)
  
    t1 <- k * (log(sf) + y + log(1 + minMu/exp(y)))
    t2 <- (k + theta) * (log(sf) + y + log(1 + minMu/exp(y))  + log(1+theta/(sf * (minMu + exp(y))))  )
    ll <- (t1 - t2)
    ll <- mean(ll[exclusionMask])
    
    # if(!is.finite(ll) & debugMyCode){
    #   browser()
    # }
  
    return(-ll)
}


lossGradENonOutlier <- function(e, D, k, b, x, sf, theta, minMu=0.00, exclusionMask, ...){
    E <-matrix(e, nrow=ncol(k))
    theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    # dW:
    y <- t(t(x %*% E %*% t(D)) + b)

    yexp <- sf * (minMu + exp(y))
    #k1 <- k *sf* exp(y) / yexp         
    k1 <- k / (1 + minMu/exp(y) ) 
    k1[exclusionMask] <- 0
    #kt <- (k + theta) *sf* exp(y) / (yexp + theta)
    kt <- (k + theta) / ( 1 + (minMu + theta/sf)/exp(y) )
    kt[exclusionMask] <- 0
    t1 <- t(x) %*% (k1 %*% D)
    t3 <- t(x) %*% (kt %*% D)
    
    # answers dE 
    dE <- (-t1 + t3)/sum(exclusionMask==FALSE)
    
    return(dE)
}
