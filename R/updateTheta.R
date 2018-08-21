updateTheta <- function(ods, thetaRange, BPPARAM, CoxR=FALSE){
    normalizationFactors(ods) <- t(predictED(ods=ods))
    mu <- normalizationFactors(ods)
    cts <- counts(ods)
    H <- x(ods) %*% E(ods)
    if('thetaCorrection' %in% names(colData(ods))){
        thetaC <- colData(ods)[['thetaCorrection']]
    } else {
        thetaC <- 1
    }
    if(isTRUE(CoxR)){
        nll <- negCRlogLikelihood
    }else{
        nll <- neglogLikelihood
    }
    
    fitparameters <- bplapply(seq_along(ods), estTheta, mu=mu, H=H,
            cts=cts, exclusionMask=NULL, thetaRange=thetaRange,
            BPPARAM=BPPARAM, nll=nll, thetaC=thetaC)
    
    theta(ods) <- vapply(fitparameters, "[[", double(1), "minimum")
    
    validObject(ods)
    
    print(summary(theta(ods)))
    return(ods)
}

thetaCorrection <- function(ods){
    mu <- normalizationFactors(ods)
    k <- counts(ods)
    
    pointVar <- (k - mu)^2
    thetaMat <- mu^2/(pointVar - mu)
    medianTheta <- rowMedians(thetaMat)
    
    thetaRatios <- thetaMat/medianTheta
    colData(ods)[['thetaCorrection']] <- colMeans(thetaRatios)
    validateOutriderDataSet(ods)
    return(ods)
}


estTheta <- function(index, cts, mu, H, thetaC, thetaRange, exclusionMask, nll){
    ctsi <- cts[index,]
    if(is.matrix(mu)){
        mu <- mu[index,]
    }
    stopifnot(!is.null(mu))
    
    # needs to be fixed in case ever used with H.
    if(!is.null(exclusionMask)){
        ctsi <- ctsi[!exclusionMask[index,]]
        mu <- mu[!exclusionMask[index,]]
    }
    est <- optimize(f= nll, interval = thetaRange, k=ctsi, mu=mu,
                    H=H, thetaC=thetaC)
} 


neglogLikelihood <- function(theta, k, mu, H, thetaC){
    theta <- theta * thetaC
    -sum(dnbinom(x=k, size=theta, mu=mu, log=TRUE))
}

negCRlogLikelihood <- function(theta, k, mu, H, thetaC){
    theta <- theta * thetaC
    H <- cbind(1,H)
    # construct a diagonal matrix.
    w <- matrix(0, ncol=nrow(H), nrow=nrow(H))
    diag(w) <- 1/(1/mu + 1/theta)
    -sum(dnbinom(x=k, size=theta, mu=mu, log=TRUE)) + 
        0.5*log(det(t(H) %*% w %*% H))
}

## grad theta is still in the fitNB file.
gradnegCRlogLikelihood <- function(theta, k, mu, H, thetaC){
    theta <- theta * thetaC
    H <- cbind(1,H)
    
    # construct a diagonal matrix.
    w <- matrix(0, ncol=nrow(H), nrow=nrow(H))
    diag(w) <- mu * theta/(mu + theta)
    w2 <- matrix(0, ncol=nrow(H), nrow=nrow(H))
    diag(w2) <- mu^2/(mu + theta)^2
    CR <- t(H) %*% w %*% H
    gradCR <- t(H) %*% w2 %*% H
    gradTheta(theta, k, mu) + 
        0.5/(det(CR)) *
        sum(diag(adjoint(CR)%*%gradCR))
}