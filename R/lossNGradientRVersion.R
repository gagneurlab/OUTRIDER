#' 
#' This file contains all r implementations of the loss functions.
#' They are only here for debug perpose and are not used within the package.
#' 
#' @noRd 
NULL

lossD <- function(par, k, H, sf, theta){
    b <- par[1]
    d <- par[-1]
    
    y <- H %*% d + b
    yexp <- sf * exp(y)
    
    ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    return(-ll)
}

lossDtrunc <- function(par, k, H, sf, theta, thetaC){
    b <- par[1]
    d <- par[-1]
    theta <- thetaC * theta
    
    y <- H %*% d + b
    yexp <- sf * exp(y)
    
    #ll = mean(k * log(yexp) - (k + theta)*log(yexp + theta))
    
    t1 <- k * (log(sf) + y) 
    t2 <- (k + theta) * (log(sf) + y + log(1+theta/(sf * exp(y))))
    ll <- mean(t1 - t2)
    
    return(-ll)
}

gradD <- function(par, k, H, sf=1, theta, thetaC){
    b <- par[1]
    d <- par[-1]
    
    theta <- thetaC * theta
    
    y <- c(H %*% d + b)
    yexp <- sf * exp(y)
    
    t1 <- colMeans(k * H)
    kt <- (k + theta) / (1 + theta/(sf*exp(y)))
    t2 <- colMeans(kt * H) 
    
    dd <- t2 - t1
    db <- mean(kt - k)
    
    return(c(db, dd))
}

lossE <- function(e, D, k, b, x, sf, theta, thetaC){
    E <- matrix(e, nrow=ncol(k))
    thetaMat <- outer(thetaC, theta)
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=t(thetaMat), log=TRUE)
    - mean( ll )
}

lossEtruncNonOutlier <- function(e, D, k, b, x, sf, theta, thetaC, 
                    exclusionMask){
    E <-matrix(e, nrow=ncol(k))
    thetaMat <- outer(thetaC, theta)
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y <- apply(y, 2, pmax, -700)
    
    t1 <- k * (log(sf) + y)
    t2 <- (k + thetaMat) * (log(sf) + y + log(1+thetaMat/(sf * exp(y))))
    
    ll <- mean((t1 - t2) * exclusionMask)
    return(-ll)
}

lossGradENonOutlier <- function(e, D, k, b, x, sf, theta, thetaC, 
                    exclusionMask){
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
