#'
#' Update D function
#' 
#' @noRd
updateDTheta <- function(ods, lasso, control, BPPARAM, optim=TRUE){
    D <- D(ods)
    b <- b(ods)
    H <- H(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    mask <- exclusionMask(ods)
    theta <- theta(ods)
    thetaC <- 1
    lambda <- lambda(ods)
    
    fitls <- bplapply(1:nrow(ods), singleDThetaFit, D=D, b=b, k=k, sf=sf, H=H, 
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
    theta(ods) <- exp(parMat[1,])
    b(ods) <- parMat[2,]
    D(ods) <- t(parMat)[,-1:-2]
    
    metadata(ods)[['Dfits']] <- fitls
    
    return(ods)
}


singleDThetaFit <- function(i, D, b, k, theta, mask, lasso, lambda, optim, control, ...){
    pari <- c(log(theta[i]), b[i], D[i,])
    ki <- k[,i]
    maski <- mask[i,]
    lambdai <- lambda[i]
    
   
    fit <- lbfgs(lossDTheta, gradDTheta,pari, k=ki, 
                ..., orthantwise_c=lambdai, 
                orthantwise_start=2,orthantwise_end = length(pari), invisible=0)
    return(fit)
}

lossDTheta <- function(par, k, H, sf, ...){
    theta <- exp(par[1])
    b <- par[2]
    d <- par[-1:-2]
    
    y <- H %*% d + b
    yexp <- sf * exp(y)
    
    ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    return(-ll)
}


gradDTheta <- function(par, k, H, sf=1, thetaC){
    theta <- exp(par[1])
    b <- par[2]
    d <- par[-1:-2]
    
    theta <- thetaC * theta
    
    y <- c(H %*% d + b)
    yexp <- sf * exp(y)
    
    t1 <- colMeans(k * H)
    kt <- (k + theta) / (1 + theta/(sf*exp(y)))
    t2 <- colMeans(kt * H) 
    
    dd <- t2 - t1
    db <- mean(kt - k)
    dtheta <- gradTheta(theta, k, yexp)/theta
    return(c(dtheta, db, dd))
}



