#'
#' Update E step for the autoencoder fit
#' 
#' @noRd
updateE <- function(ods, control, BPPARAM, thetaCorrection=FALSE){
    e <- as.vector(E(ods))
    D <- D(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    x <- x(ods)
    b <- b(ods)
    theta <- theta(ods)
    mask <- t(exclusionMask(ods))
    if(isTRUE(thetaCorrection)){
        thetaC <- colData(ods)[['thetaCorrection']]
    } else {
        thetaC <- rep(1,ncol(ods))
    }
    
    control$trace <- 3
    fit <- optim(e, fn=truncLogLiklihoodE, gr=gradientE,
            k=k, x=x, sf=sf, D=D, b=b, theta=theta, exclusionMask=mask, 
            thetaC=thetaC,
            method="L-BFGS-B", lower=-100, upper=100, control=control)
    
    # Check that fit converged
    if(fit$convergence!=0){
        print(paste('Update E did not converge: ', fit$message))
    }
    
    # update ods
    E(ods) <- fit$par
    
    return(ods)
}

lossE <- function(e, D, k, b, x, sf, theta, thetaC){
    
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    theta <- outer(thetaC, theta)
    ## encoding 
    E <-matrix(e, nrow=ncol(k))
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=t(theta), log=TRUE)
    - mean( ll )
}

lossEtrunc <- function(e, D, k, b, x, sf, theta, thetaC){
    theta <- outer(thetaC, theta)
    E <-matrix(e, nrow=ncol(k))
    
    y <- t(t(x %*% E %*% t(D)) + b)
    
    #ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    #ll = mean(k * log(yexp) - (k + theta)*log(yexp + theta))
    
    t1 <- k * (log(sf) + y)
    t2 <- (k + theta) * (log(sf) + y + log(1+theta/(sf * exp(y)))  )
    ll <- mean(t1 - t2)
    
    # if(!is.finite(ll) & debugMyCode){
    #   browser()
    # }
    
    return(-ll)
}


lossEtruncNonOutlier <- function(e, D, k, b, x, sf, theta, thetaC, exlusionMask){
    theta <- outer(thetaC, theta)
    E <-matrix(e, nrow=ncol(k))
  
    y <- t(t(x %*% E %*% t(D)) + b)
  
    t1 <- k * (log(sf) + y)
    t2 <- (k + theta) * (log(sf) + y + log(1+theta/(sf * exp(y)))  )
    ll <- (t1 - t2)
    ll <- mean(ll[exclusionMask])
    
    # if(!is.finite(ll) & debugMyCode){
    #   browser()
    # }
  
    return(-ll)
}


lossGradENonOutlier <- function(e, D, k, b, x, sf, theta, thetaC, exclusionMask){
    theta <- outer(thetaC, theta)
    E <-matrix(e, nrow=ncol(k))
    theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    # dW:
    y <- t(t(x %*% E %*% t(D)) + b)

    yexp <- sf * exp(y)
    #k1 <- k *sf* exp(y) / yexp         
    kt <- (k + theta) / ( 1 + theta/(sf*exp(y)) )
    k[exclusionMask] <- 0
    kt[exclusionMask] <- 0
    t1 <- t(x) %*% (k %*% D)
    t3 <- t(x) %*% (kt %*% D)
    
    # answers dE 
    dE <- (-t1 + t3)/sum(exclusionMask==FALSE)
    
    return(dE)
}
