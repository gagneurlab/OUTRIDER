#'
#' Update D function
#' 
#' @noRd
updateD <- function(ods, lasso, control, BPPARAM, optim=TRUE){
    D <- D(ods)
    b <- b(ods)
    H <- H(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    mask <- exclusionMask(ods)
    theta <- theta(ods)
    thetaC <- thetaCorrection(ods)
    lambda <- lambda(ods)
    
    fitls <- bplapply(1:nrow(ods), singleDFit, D=D, b=b, k=k, sf=sf, H=H, 
            theta=theta, mask=mask, control=control, thetaC=thetaC, 
            lasso=lasso, lambda=lambda, BPPARAM=BPPARAM, optim=optim)
    
    # update D and bias terms
    parMat <- sapply(fitls, '[[', 'par')
    if(isTRUE(optim)){
        mcols(ods)['FitDMessage'] <- sapply(fitls, '[[', 'message')
        print(table(mcols(ods)[,'FitDMessage']))
        mcols(ods)[,'NumConvergedD'] <- mcols(ods)[,'NumConvergedD'] + grepl(
        "CONVERGENCE: REL_REDUCTION_OF_F .. FACTR.EPSMCH", 
        mcols(ods)[,'FitDMessage'])
    }else{
        ## eventually add Num converged for OWL-QN...
    }
    b(ods) <- parMat[1,]
    D(ods) <- t(parMat)[,-1]
    
    metadata(ods)[['Dfits']] <- fitls
    
    return(ods)
}


singleDFit <- function(i, D, b, k, theta, mask, lasso, lambda, optim, control, ...){
    pari <- c(b[i], D[i,])
    ki <- k[,i]
    thetai <- theta[i]
    maski <- mask[i,]
    lambdai <- lambda[i]
    
    if(isTRUE(lasso)){
        if(isTRUE(optim)){
            fit <- optim(pari, fn=truncLogLiklihoodDLasso, gr=gradientDLasso, 
                     k=ki, theta=thetai, exclusionMask=maski, lambda=lambdai,
                     ..., lower=-100, upper=100, method='L-BFGS', control=control)
        } else {
        fit <- lbfgs(truncLogLiklihoodD, gradientD,pari, k=ki, theta=thetai, 
                     exclusionMask=maski, ..., orthantwise_c=lambdai, 
                     orthantwise_start=1,orthantwise_end = length(pari), invisible=1)
        }
        
    }else{
        fit <- optim(pari, fn=truncLogLiklihoodD, gr=gradientD, 
                     k=ki, theta=thetai, exclusionMask=maski, ..., control=control,
                     lower=-100, upper=100, method='L-BFGS')   
    }
    return(fit)
}

lossD <- function(par, k, H, sf, theta){
    b <- par[1]
    d <- par[-1]
    
    y <- H %*% d + b
    yexp <- sf * exp(y)
    
    ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    return(-ll)
}

predictGeneMu <- function(par, k, H, sf){
    b <- par[1]
    d <- par[-1]
    
    y <- H %*% d + b
    yexp <- sf * exp(y)
    return(yexp)
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
