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
    thetaC <- thetaCorrection(ods)
    
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
    E <-matrix(e, nrow=ncol(k))
    thetaMat <- outer(thetaC, theta)
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=t(thetaMat), log=TRUE)
    - mean( ll )
}

lossEtruncNonOutlier <- function(e, D, k, b, x, sf, theta, thetaC, exclusionMask){
    E <-matrix(e, nrow=ncol(k))
    thetaMat <- outer(thetaC, theta)
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y <- apply(y, 2, pmax, -700)
    
    t1 <- k * (log(sf) + y)
    t2 <- (k + thetaMat) * (log(sf) + y + log(1+thetaMat/(sf * exp(y))))
    
    ll <- mean((t1 - t2) * exclusionMask)
    return(-ll)
}

lossGradENonOutlier <- function(e, D, k, b, x, sf, theta, thetaC, exclusionMask){
    E <-matrix(e, nrow=ncol(k))
    thetaMat <- outer(thetaC, theta)
    
    # dW:
    y <- t(t(x %*% E %*% t(D)) + b)
    y <- apply(y, 2, pmax, -700)

    kt <- (k + thetaMat) / ( 1 + thetaMat/(sf*exp(y)) )
    
    t1 <- t(x) %*% ((k  * exclusionMask) %*% D)
    t3 <- t(x) %*% ((kt * exclusionMask) %*% D)
    
    # answers dE 
    dE <- (-t1 + t3)/sum(exclusionMask)
    
    return(dE)
}

if(FALSE){
    fn <- truncLogLiklihoodE
    gr <- gradientE
    fn <- lossEtruncNonOutlier
    gr <- lossGradENonOutlier
    fit <- optim(e, fn=fn, gr=gr,
                 k=k, x=x, sf=sf, D=D, b=b, theta=theta, exclusionMask=mask, 
                 thetaC=thetaC,
                 method="L-BFGS-B", lower=-100, upper=100, control=control)
    
}
