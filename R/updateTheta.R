
#' 
#' Update theta step for autoencoder
#' 
#' @noRd
updateTheta <- function(ods, thetaRange, BPPARAM){
    normalizationFactors(ods) <- t(predictC(ods))
    mu <- normalizationFactors(ods)
    cts <- counts(ods)
    H <- H(ods)
    thetaC <- thetaCorrection(ods)
    sMask <- sampleExclusionMask(ods, aeMatrix=TRUE)
    
    fitparameters <- bplapply(seq_along(ods), estTheta, mu=mu, H=H,
            cts=cts, exclusionMask=sMask, thetaRange=thetaRange,
            BPPARAM=BPPARAM, nll=negLogLikelihoodTheta, thetaC=thetaC)
    
    theta(ods) <- vapply(fitparameters, "[[", double(1), "minimum")
    print(summary(theta(ods)))
    
    validObject(ods)
    return(ods)
}

estTheta <- function(index, cts, mu, H, thetaC, thetaRange, exclusionMask, nll){
    sMaski  <- exclusionMask[index,]
    ctsi    <- cts[index,sMaski]
    mui     <- mu[index,sMaski]
    thetaCi <- thetaC[sMaski]
    
    est <- optimize(f=nll, 
            interval=thetaRange, k=ctsi, mu=mui, thetaC=thetaCi)
} 

negLogLikelihoodTheta <- function(theta, k, mu, thetaC){
    theta <- theta * thetaC
    -sum(dnbinom(x=k, size=theta, mu=mu, log=TRUE))
}
