#'
#' Main autoencoder fit function
#'
fitAutoencoder <- function(ods, q, thetaRange=c(1e-2, 1e3), 
                    convergence=1e-5, iterations=15, initialize=TRUE,
                    correctTheta='none', usePCA=TRUE, lasso=FALSE, 
                    runLassoFit=TRUE, nFolds=20,
                    control=list(), useOptim=TRUE, L1encoder=FALSE,
                    newCVversion=FALSE, useSE=FALSE,
                    BPPARAM=bpparam()){
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
        ods <- initAutoencoder(ods, q, thetaRange, usePCA=usePCA)
    }
     
    # initial loss
    lossList <- c(init_pca=lossED(ods))
    print(paste0('Initial PCA loss: ', lossList[1]))
    convList <- numeric()
    
    #' initialize D 
    ods <- updateD(ods, lasso=lasso, control=control, BPPARAM=BPPARAM, 
            optim=useOptim)
    lossList <- updateLossList(ods, lossList, 'init', 'D')
    
    # initialize theta step
    ods <- updateTheta(ods, thetaRange, correctTheta=correctTheta, BPPARAM=BPPARAM)
    lossList <- updateLossList(ods, lossList, 'init', 'Theta')
    
    # optimize log likelihood
    t1 <- Sys.time()
    currentLoss <- lossED(ods)
    for(i in seq_len(iterations)){
        t2 <- Sys.time()
        
        
        # update lasso
        if(isTRUE(lasso) & i == 2 & isTRUE(runLassoFit)){
            ods <- updateLambda(ods, nFolds=nFolds, control=control, 
                                BPPARAM=BPPARAM, optim=useOptim, newCVversion=newCVversion,
                                useSE=useSE)
        }
        
        # update E step
        ods <- updateE(ods, control=control, BPPARAM=BPPARAM, L1encoder=L1encoder)
        lossList <- updateLossList(ods, lossList, i, 'E')
        convList <- updateConvergenceList(convList, lossList, ods, L1encoder=L1encoder)
        
        # update D step
        ods <- updateD(ods, lasso=lasso, control=control, BPPARAM=BPPARAM, optim=useOptim)
        lossList <- updateLossList(ods, lossList, i, 'D')
        convList <- updateConvergenceList(convList, lossList, ods, L1encoder=L1encoder)
        
        # update theta step
        ods <- updateTheta(ods, thetaRange, correctTheta=correctTheta, BPPARAM=BPPARAM)
        lossList <- updateLossList(ods, lossList, i, 'theta')
        convList <- updateConvergenceList(convList, lossList, ods, L1encoder=L1encoder)
        
        print(paste('Time for one autoencoder loop:', Sys.time() - t2))
        
        # check 
        if(all(abs(currentLoss - convList[length(convList) - 2:0]) < convergence)){
            message(date(), ': the AE correction converged with:',
                    lossList[length(lossList)])
            break
        }
        currentLoss <- convList[length(convList)]
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
    metadata(ods)[['convList']] <- convList
    
    validObject(ods)
    return(ods)
}

initAutoencoder <- function(ods, q, thetaRange, usePCA){
    
    pca <- pca(x(ods), nPcs=q)
    pc  <- loadings(pca)
    
    if(isFALSE(usePCA)){
        rnorm(length(pc), sd=sd(pc))
    }
    
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

updateConvergenceList <- function(convList, lossList, ods, L1encoder){
    currLoss <- lossList[length(lossList)] + sum(lambda(ods)*rowSums(abs(D(ods))))
    if(isTRUE(L1encoder)){
        currLoss <- currLoss + mean(lambda(ods)) * sum(abs(E(ods)))
    }
    convList <- c(convList, currLoss)
    return(convList)
}

lossED <- function(ods, step=c('none', 'E', 'D', 'Theta')){
    k <- t(counts(ods))
    y_exp <- predictC(ods)
    theta <- outer(thetaCorrection(ods), theta(ods))
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=t(theta), log=TRUE)
    
    # compute lasso regularization
    w <- switch(match.arg(step),
                none = { 0 },
                E = { sum(mean(lambda(ods)) * abs(E(ods))) },
                D = { sum(lambda(ods) * abs(D(ods))) },
                Theta = { 0 },
                stop('Please provide a correct step for lossED.')
    )
    
    return( - (mean(ll) + w) )
}
