#'
#' Update E step for the autoencoder fit
#' 
#' @noRd
updateE <- function(ods, control, BPPARAM, verbose){
    e <- as.vector(E(ods))
    D <- D(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    x <- x(ods)
    b <- b(ods)
    theta <- theta(ods)
    sMask <- t(sampleExclusionMask(ods, aeMatrix=TRUE))
    thetaC <- thetaCorrection(ods)
    
    fit <- optim(e, fn=truncLogLiklihoodE, gr=gradientE, k=k, x=x,
            sf=sf, D=D, b=b, theta=theta, exclusionMask=sMask, thetaC=thetaC,
            method="L-BFGS-B", lower=-100, upper=100, control=control)
    
    # Check that fit converged
    if(isTRUE(verbose) & fit$convergence != 0){
        print(paste('Update E did not converge: ', fit$message))
    }
    
    # update ods
    E(ods) <- fit$par
    
    return(ods)
}
