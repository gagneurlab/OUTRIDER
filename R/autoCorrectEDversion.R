#'
#' Main autoencoder fit function
#'
fitAutoencoder <- function(ods, q, robust=TRUE, thetaRange=c(0.1, 250), 
                    convergence=1e-5, loops=15, pValCutoff=0.0,
                    initialize=TRUE, noRobustLast=TRUE, coxReid=FALSE, 
                    thetaCorrection=FALSE,
                    control=list(), BPPARAM=bpparam(), ...){
    
    # Check input
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    checkSizeFactors(ods)
    if(!bpisup(BPPARAM)){
        bpstart(BPPARAM)
    }
    # reset counters 
    mcols(ods)['NumConvergedD'] <- 0
    
    k <- t(counts(ods, normalized=FALSE))
    sf <- sizeFactors(ods)
    
    # initialize W using PCA and bias as zeros.
    if(isTRUE(initialize) | is.null(E(ods)) | is.null(D(ods))){
        ods <- initAutoencoder(ods, q, thetaRange)
    }
     
    # initial loss
    lossList <- lossED(ods)
    print(paste0('Initial PCA loss: ', lossList[1]))
    
    
    # optimize log likelihood
    t1 <- Sys.time()
    currentLoss <- lossED(ods)
    for(i in seq_len(loops)){
        t2 <- Sys.time()
        
        # update D step
        ods <- updateD(ods, control=control, BPPARAM=BPPARAM)
        lossList[i*3-1] <- lossED(ods)
        print(paste0('Iteration: ', i, '; update D loss: ', lossList[i*3-1]))
        
        # update Theta correction
        
        # update theta step
        ods <- updateTheta(ods, thetaRange, CoxR=coxReid, BPPARAM)
        lossList[i*3+0] <- lossED(ods)
        print(paste0('Iteration: ', i, ' theta loss: ', lossList[i*3+0]))
        
        
        # mask outlier
        if(i != 1 & isTRUE(robust)){
            ods <- maskOutliers(ods, pValCutoff=pValCutoff, BPPARAM=BPPARAM)
        } else {
            exclusionMask(ods) <- 1
        }
        
        # update E step
        ods <- updateE(ods, control=control, BPPARAM=BPPARAM)
        lossList[i*3+1] <- lossED(ods)
        print(paste0('Iteration: ', i, ' update E loss: ', lossList[i*3+1]))
        
        
        print(paste('Time for one autoencoder loop:', Sys.time() - t2))
        
        # check 
        if(all(abs(currentLoss - lossList[i*3+c(-1,0,1)]) < convergence)){
            message(date(), ': the AE correction converged with:',
                    lossList[i*3+1])
            break
        }
        currentLoss <- lossList[i*3+1]
    }
    
    # Final model fit
    if(isTRUE(noRobustLast)){
        exclusionMask(ods) <- 1
    }
    ods <- updateD(ods, control, BPPARAM)
    # no Inf with optimize.
    ods <- updateTheta(ods, c(0, 500), BPPARAM)
    
    print(Sys.time() - t1)
    print(paste0(i, ' Final nb-AE loss: ', lossED(ods)))
    
    bpstop(BPPARAM)
    
    # add correction factors
    correctionFactors <- t(predictED(ods=ods))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    normalizationFactors(ods) <- correctionFactors
    
    # add additional values for the user to the object
    metadata(ods)[['dim']] <- dim(ods)
    metadata(ods)[['loss']] <- lossList
    
    validObject(ods)
    return(ods)
}

initAutoencoder <- function(ods, q, thetaRange, BPPARAM){
    pca <- pca(x(ods), nPcs=q)
    pc  <- loadings(pca)
    
    # Set initial values from PCA
    D(ods) <- pc
    E(ods) <- pc
    b(ods) <- rowMeans(log(counts(ods) + 1))
    
    # initialize theta
    theta(ods) <- robustMethodOfMomentsOfTheta(counts(ods), 
            minTheta=thetaRange[1], maxTheta=thetaRange[2])
    
    return(ods)
}


maskOutliers <- function(ods, pValCutoff=0.01, BPPARAM){
    ods <- computePvalues(ods, BPPARAM=BPPARAM, method='None')
    mask <- matrix(1, nrow=nrow(ods), ncol=ncol(ods))
    mask[pValue(ods) < pValCutoff/ncol(ods)] <- 0
    
    print(paste(sum(mask==0), 'outliers identified in this iteration.'))
    print(paste('Top 5 masked genes: ', paste(collapse=", ",
            tail(sort(rowSums(mask == 0))))))
    
    exclusionMask(ods) <- mask
    return(ods)
}

lossED <- function(ods){
    
    ## encoding 
    b <- b(ods)
    E <- E(ods)
    D <- D(ods)
    x <- x(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    theta <- theta(ods)
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=theta, log=TRUE)
    - mean( ll )
}

predictED <- function(ods){
    E <- E(ods)
    D <- D(ods)
    b <- b(ods)
    x <- x(ods)
    sf <- sizeFactors(ods)

    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * exp(y)
    y_exp
}
