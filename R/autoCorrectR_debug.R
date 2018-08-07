robThetaFade200L20It25 <- function(ods, q, debug=FALSE, ...){
    autoCorrectRCooksIter2Debug(ods, q=q,
            robust='iterative', useDESeq=FALSE,
            modelTheta='fade', initTheta=200, trim=0.05,
            mask=TRUE, loops=20, internIter=25, ...)
}

robMix25L5I40 <- function(ods, q, debug=FALSE, ...){
    autoCorrectRCooksIter2Debug(ods, q=q, robust='iterative', useDESeq=FALSE,
            modelTheta='mix', initTheta=25, mask=TRUE, loops=5, internIter=40, ...)
}


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

autoCorrectRCooksIter2Debug <- function(ods, q, initTheta=25, maskOnce=FALSE,
                    robust=c('once', 'iterative', 'none'), pcaOnly=FALSE,
                    modelTheta=c('no', 'fit', 'mean', 'fade', 'fadeM5', 'trimmed', 't2_25', 't1_10', 't0_5', 'tw0_5', 'mix', 'fix'),
                    internIter=10, loops=10, debug=TRUE, trim=0, useDESeq=TRUE,
                    mask=FALSE, control=list(), cLoss = TRUE, BPPARAM=bpparam(),
                    ThetaCooks=FALSE, pValTest=FALSE, poisson='no', ...){
    methodStr <- paste(
        'q:', q, 'initTheta:', round(mean(initTheta), 2), 'robust:', robust,
        'pcaOnly:', pcaOnly, 'modelTheta:', modelTheta, 'useDESeq:', useDESeq,
        'internIter:', internIter, 'loops:', loops, 'trim:', trim, 
        'cLoss:', cLoss, 'ThetaCooks:', ThetaCooks)
    message('Using: ', methodStr)
    
    # set defaults
    robust <- match.arg(robust)
    modelTheta <- match.arg(modelTheta)
    curMask <- FALSE
    theta <- initTheta
    myLoss <- loss2
    myLossGrad <- lossGrad2
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    if(isScalarNumeric(theta)){
        theta <- rep(initTheta, nrow(ods))
    }
    if(isTRUE(mask)){
        myLoss <- function(w, k, x, s, xbar, theta, outlier, poisson='no') { lossNonOutlier(w, k, x, s, xbar, theta, outlier) }
        myLossGrad <- function(w, k, x, s, xbar, theta, outlier, poisson='no') { lossGradNonOutlier(w, k, x, s, xbar, theta, outlier) }
        if(isTRUE(cLoss)){
            myLoss <- function(w, k, x, s, xbar, theta, outlier, poisson='no') { lossNonOutlierC(w, k, x, s, xbar, theta, outlier) }
            myLossGrad <- function(w, k, x, s, xbar, theta, outlier, poisson='no') { lossGradNonOutlierC(w, k, x, s, xbar, theta, outlier) }
        }
    }
    if(poisson != 'no'){
        myLoss <- loss2
        myLossGrad <- lossGrad2
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    k_no <- k
    rep_k <- NULL
    
    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
    xbar <- apply(x0, 2, mean, trim=trim)
    x <- t(t(x0) - xbar)
    
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    
    # check initial loss
    print(paste0(date(), ': Initial PCA loss: ',
            myLoss(w_guess, k, x, s, xbar, theta, poisson=poisson, 
                    outlier=curMask))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    
    w_fit <- w_guess
    for(i in seq_len(loops)){
        if(isTRUE(pcaOnly)){
            fit <- list(par=w_fit)
            break
        }
        
        if(robust == 'iterative' || robust == 'once' & i == 1){
            if(isTRUE(pValTest)){
                mu <- predictC(w_fit, k, s, xbar)
                rep_k <- findOutlierNBfit(k, mu, pValCutoff=0.01)
                k_no <- k
            }else{
                mu <- predictC(w_fit, k, s, xbar)
                rep_k <-replaceOutliersCooks(k, mu, q=q,
                        BPPARAM=BPPARAM, theta=modelTheta != 'no', useDESeq=useDESeq,
                        ThetaCooks=ThetaCooks)
            
                k_no <- rep_k[['cts']]
            }
            if(modelTheta != 'no'){
                theta <- getModeledTheta(modelTheta, rep_k[['theta']], theta, k=k, mu=mu, i, loops)
                message(paste(round(summary(theta), 2), names(summary(theta)), collapse=', '))
            }
            if(isTRUE(mask)){
                k_no <- k
                if(isFALSE(maskOnce) | i == 1){
                    curMask <- rep_k$mask
                } else {
                    curMask <- curMask | rep_k$mask
                    message('Merge masked: ', sum(curMask))
                }
                assay(ods, 'excludeMask') <- t(curMask)
            }
        }
        
        x0 <- log((1+k_no)/s)
        x <- t(t(x0) - xbar)
        
        control$maxit <- internIter
        fit <- autoCorrectFit(w_fit, loss=myLoss, lossGrad=myLossGrad, 
                poisson=poisson, 
                k_no, x, s, xbar, theta, outlier=curMask, control, ...)
        
        w_fit <- fit$par
        message(date(), ': Iteration ', i, ' loss: ', 
            signif(myLoss(w_fit,   k_no, x, s, xbar, theta, poisson=poisson, outlier=curMask), 8),
            ' PCA-loss: ', 
            signif(myLoss(w_guess, k_no, x, s, xbar, theta, poisson=poisson, outlier=curMask), 8))
        
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
        
        if(isTRUE(debug)){
            metadata(ods)[[paste0('iter_', i)]] <- list(
                w = w_fit,
                loss = myLoss(w_fit, k_no, x, s, xbar, theta, poisson=poisson, outlier=curMask),
                lossGrad = myLossGrad(w_fit, k_no, x, s, xbar, theta, poisson=poisson, outlier=curMask),
                fit=fit,
                dds=repdds,
                ods=repods
            )
        }
    }
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               myLoss(w_fit,k, x, s,xbar, theta, poisson=poisson, outlier=curMask))
    )
    
    correctionFactors <- t(predictC(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    
    message('Used: ', methodStr)
    return(ods)
}   

getModeledTheta <- function(modelTheta=c('fit', 'mean', 'fade', 'fadeM5', 'trimmed', 't2_25', 't1_10', 't0_5', 'tw0_5', 'mix', 'fix'), 
                    fitT, oldT, k, mu, it, loops){
    ans <- switch(match.arg(modelTheta),
        fit     = { fitT },
        mean    = { (fitT + oldT)/2 },
        fade    = { (mean(oldT)*(loops - it + 1) + fitT*(it))/(loops + 1) },
        fadeM5  = { pmin(1000, pmax(5, (mean(oldT)*(loops - it + 1) + fitT*(it))/(loops + 1))) },
        trimmed = { pmin(pmax((fitT + oldT) / 2, 2), 1000) },
        t2_25   = { pmin(25, pmax(2, (mean(oldT)*(loops - it + 1) + fitT*(it))/(loops + 1))) },
        t1_10   = { pmin(10, pmax(1, (mean(oldT)*(loops - it + 1) + fitT*(it))/(loops + 1))) },
        t0_5    = { pmin(5, pmax(0.1, (oldT + fitT)/2)) },
        tw0_5   = { pmin(5, pmax(0.5, (oldT + fitT*2)/3)) },
        mix     = { pmin(5000, pmax(1.5, (mean(oldT)*(loops - it + 1)/2 + fitT*(it + 1)*2)/((loops + 3*it + 5)/2)))},
        fix     = { pmin(50000, pmax(1.1, fitT + 0.05*(fitT)^2)) },
        stop('Option not known: ', modelTheta)
    )
    return(ans)
}
