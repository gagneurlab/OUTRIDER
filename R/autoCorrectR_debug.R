


getExpectations <- function(w, k, s, xbar){
    t(predictC(w, k, s, xbar))
}

getAEData <- function(ods, w, replace=FALSE, BPPARAM=SerialParam()){
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    if(isTRUE(replace)){
        k <-replaceOutliersCooks(k, BPPARAM=BPPARAM)
    }
    
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    ans <- list(k=k, s=s, xbar=xbar)
    if(!missing(w)){
        ans["norm"] <- list(t(predictC(w, k, s, xbar)))
    }
    
    return(ans)
}

autoCorrectRCooksIter2Debug <- function(ods, q, initTheta=25, modelTheta=FALSE,
                    robust=c('once', 'iterative', 'none'), pcaOnly=FALSE,
                    internIter=10, noFirst=FALSE,
                    control=list(), BPPARAM=bpparam(), ...){
    robust <- match.arg(robust)
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    myLoss <- loss
    myLossGrad <- lossGrad
    theta <- initTheta
    if(isTRUE(modelTheta)){
        myLoss <- loss2
        myLossGrad <- lossGrad2
        theta <- rep(initTheta, nrow(ods))
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    rep_k <- NULL
    k_no <- k
    
    if(robust != 'none' & isFALSE(noFirst)){
        k_no <-replaceOutliersCooks(k_no, BPPARAM=BPPARAM)
        if(isTRUE(modelTheta)){
            theta <- k_no[['theta']]
            k_no  <- k_no[['cts']]
        } else {
            k_no  <- k_no[['cts']]
        }
    }
    
    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0(date(), ': Initial PCA loss: ',
               myLoss(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    
    w_fit <- w_guess
    for(i in 1:10){
        if(isTRUE(pcaOnly)){
            fit <- list(par=w_fit)
            break
        }
        
        if(robust == 'iterative' || robust == 'once' & i == 1){
            rep_k <-replaceOutliersCooks(k,predictC(w_fit, k, s, xbar), 
                    BPPARAM=BPPARAM, theta=modelTheta)
            if(isTRUE(modelTheta)){
                theta <- rep_k[['theta']]
                k_no <- rep_k[['cts']]
            } else {
                k_no <- rep_k[['cts']]
            }
        } else if(robust == 'once'){
            k_no <- k_no
        } else {
            k_no <- k
        }
        
        x0 <- log((1+k_no)/s)
        x <- t(t(x0) - xbar)
        
        control$maxit <- internIter
        fit <- autoCorrectFit(w_fit, loss=myLoss, lossGrad=myLossGrad, k_no, x, s, xbar, theta, 
                              control, ...)
        
        w_fit <- fit$par
        message(date(), ': Iteration ', i, ' loss: ', myLoss(w_fit, k_no, x, s, xbar, theta))
        
        #Check that fit converged
        if(fit$convergence!=0){
            warning(paste0("Fit didn't converge with warning: ", fit$message))
        }
        
        repdds <- repods <- NULL 
        if('dds' %in% names(rep_k)){
            repdds <- rep_k$dds
        }
        if('ods' %in% names(rep_k)){
            repods <- rep_k$ods
        }
        
        metadata(ods)[[paste0('iter_', i)]] <- list(
            w = w_fit,
            loss = myLoss(w_fit, k_no, x, s, xbar, theta),
            lossGrad = myLossGrad(w_fit, k_no, x, s, xbar, theta),
            fit=fit,
            dds=repdds,
            ods=repods
        )
    }
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               myLoss(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}   
