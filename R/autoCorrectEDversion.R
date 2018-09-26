#' 
#' Main autoencoder fit function
#' 
#' @noRd
fitAutoencoder <- function(ods, q, thetaRange=c(1e-2, 1e3), 
                    convergence=1e-5, iterations=15, initialize=TRUE,
                    control=list(), BPPARAM=bpparam()){
    # Check input
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    checkSizeFactors(ods)
    checkThetaRange(thetaRange)
    
    if(!bpisup(BPPARAM)){
        bpstart(BPPARAM)
    }
    
    # initialize W using PCA and bias as zeros.
    if(isTRUE(initialize) | is.null(E(ods)) | is.null(D(ods))){
        ods <- initAutoencoder(ods, q, thetaRange)
    }
     
    # initial loss
    lossList <- c(init_pca=lossED(ods))
    print(paste0('Initial PCA loss: ', lossList[1]))
    
    # initialize D 
    ods <- updateD(ods, control=control, BPPARAM=BPPARAM)
    lossList <- updateLossList(ods, lossList, 'init', 'D')
    
    # initialize theta step
    ods <- updateTheta(ods, thetaRange, BPPARAM=BPPARAM)
    lossList <- updateLossList(ods, lossList, 'init', 'Theta')
    
    # optimize log likelihood
    t1 <- Sys.time()
    currentLoss <- lossED(ods)
    for(i in seq_len(iterations)){
        t2 <- Sys.time()
        
        # update E step
        ods <- updateE(ods, control=control, BPPARAM=BPPARAM)
        lossList <- updateLossList(ods, lossList, i, 'E')
        
        # update D step
        ods <- updateD(ods, control=control, BPPARAM=BPPARAM)
        lossList <- updateLossList(ods, lossList, i, 'D')
    
        # update theta step
        ods <- updateTheta(ods, thetaRange, BPPARAM=BPPARAM)
        lossList <- updateLossList(ods, lossList, i, 'theta')
        
        print(paste('Time for one autoencoder loop:', Sys.time() - t2))
        
        # check 
        curLossDiff <- abs(currentLoss - lossList[length(lossList) - 2:0])
        if(all(curLossDiff < convergence)){
            message(date(), ': the AE correction converged with:',
                    lossList[length(lossList)])
            break
        }
        currentLoss <- lossList[length(lossList)]
    }
    
    bpstop(BPPARAM)
    print(Sys.time() - t1)
    
    print(paste0(i, ' Final nb-AE loss: ', lossList[length(lossList)]))
    
    # add correction factors
    correctionFactors <- t(predictC(ods))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    normalizationFactors(ods) <- correctionFactors
    
    # add additional values for the user to the object
    metadata(ods)[['dim']] <- dim(ods)
    metadata(ods)[['loss']] <- lossList
    metadata(ods)[['convList']] <- lossList
    
    validObject(ods)
    return(ods)
}

initAutoencoder <- function(ods, q, thetaRange){
    
    pca <- pca(x(ods), nPcs=q)
    pc  <- loadings(pca)
    
    # Set initial values from PCA
    D(ods) <- pc
    E(ods) <- pc
    b(ods) <- rowMeans(log(counts(ods) + 1))
    
    # initialize theta
    theta(ods) <- robustMethodOfMomentsOfTheta(counts(ods), 
            minTheta=thetaRange[1], maxTheta=thetaRange[2])
    thetaCorrection(ods) <- 1
    
    # reset counters 
    mcols(ods)['NumConvergedD'] <- 0
    
    return(ods)
}

updateLossList <- function(ods, lossList, i, stepText){
    currLoss <- lossED(ods)
    lossList <- c(lossList, currLoss)
    names(lossList)[length(lossList)] <- paste0(i, '_', stepText)
    print(paste0(date(), ': Iteration: ', i, ' ', stepText, ' loss: ', currLoss))
    return(lossList)
}

lossED <- function(ods, step=c('none', 'E', 'D', 'Theta')){
    k <- t(counts(ods))
    y_exp <- predictC(ods)
    theta <- outer(thetaCorrection(ods), theta(ods))
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=t(theta), log=TRUE)
    
    return( - mean(ll) )
}
