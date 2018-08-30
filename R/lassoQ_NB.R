
#'
#' Lasso CV to identify lambda.
#' 
#' @noRd
updateLambda <- function(ods, nFolds=20, control, BPPARAM, optim=TRUE){
    D <- D(ods)
    b <- b(ods)
    H <- H(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    theta <- theta(ods)
    lambda = c(0, exp(c(-20, -15, seq(-10, 10, length.out=25))))
    if(isTRUE(optim)){
        lambda = c(0, exp(c(-20, -15, seq(-10, 10, length.out=25))))
    } else {
        lambda = c(0, exp(seq(-15, 0, length.out=25)))
    }
    folds <- chunk(sample(ncol(ods)), n.chunks=nFolds)
    
    t1 <- Sys.time()
    fitls <- bplapply(seq_along(ods), lassoGeneCV, D=D, b=b, k=k, sf=sf, H=H, theta=theta, 
            lambda=lambda, folds=folds, control=control, optim=optim, BPPARAM=BPPARAM)
    t2 <- Sys.time()
    
    lambda(ods) <- sapply(fitls, '[[', 'lambda')
    metadata(ods)[['lambdaFits']] <- fitls
    
    print(paste('Lasso lambda fit time: ', t2 - t1))
    print(summary(log10(lambda(ods) + exp(lambda[2])/100)))
    
    return(ods)
}


lassoGeneCV <- function(i, D, b, k, theta, lambda, H , sf, folds, optim, setPar=TRUE, debug=FALSE, ...){
    pari <- c(b[i], D[i,])
    ki <- k[,i]
    thetai <- theta[i]
    
    res <- matrix(0, nrow=length(folds), ncol=length(lambda))
    for(f in seq_along(folds)){
        pari <- c(b[i], D[i,])
        for(l in seq_along(lambda)){
            run <- lassoSingleCVrun(l=lambda[l], testSet=folds[[f]], ki=ki, H=H,
                    thetai=thetai, pari=pari, sf=sf, optim=optim)
            if(isTRUE(setPar)){
                pari <- run$par
            }
            res[f, l] <- run$loss
            #print(paste(f, round(lambda[l], 4), paste(round(sort(abs(run$par))[c(1:7, 38:45)], 4), collapse=" ")))
        }
    }
    if(isTRUE(debug)){
        par(mfrow=c(1,2))
        boxplot(lapply(1:ncol(res), function(x) pmin(1000, res[,x])), log='y')
        plot(apply(res, 2, mean, trim=0.1))
    }
    
    cmean <- round(apply(res, 2, mean, trim=0), 3)
    minIdx <- max(which(cmean == min(cmean)))
    csd <- colSds(res)
    
    return(list(
        means=cmean, 
        csd=csd, 
        minIdx=minIdx, 
        lambda=lambda[minIdx], 
        pari=pari,
        res=res))
}

lassoSingleCVrun <- function(l, testSet, ki, H, thetai, pari, sf, optim){
    if(isTRUE(optim)){
        fitpar <- optim(pari, fn=truncLogLiklihoodDLasso, gr=gradientDLasso, 
                        k=ki[-testSet], H=H[-testSet,], theta=thetai, exclusionMask=1, lambda=l,
                        sf=sf[-testSet], thetaC=1, lower=-100, upper=100, method='L-BFGS')
    } else{
        fitpar <- lbfgs(truncLogLiklihoodD, gradientD, pari,  k=ki[-testSet], H=H[-testSet,], theta=thetai, 
                        exclusionMask=1, sf=sf[-testSet], thetaC=1, orthantwise_c=l, epsilon = 1e-9,
                        orthantwise_start=1,orthantwise_end = length(pari), invisible=1)
    }
    fitpar$par[fitpar$par == -100] <- 0
    fitpar$par[fitpar$par == 100] <- 0
    list('par'=c(fitpar$par), 'loss'=lossD(fitpar$par, ki[testSet], H[testSet,], sf[testSet], thetai))
}

if(FALSE){
    ods <- odsg
    BPPARAM <- MulticoreParam(30, 150, progressbar=TRUE)
    optim <- FALSE
    nFolds <- 10
    control <- list()
    ods1 <- updateLambda(odsg, 10, list(), BPPARAM, optim=TRUE)
    ods2 <- updateLambda(odsg, 10, list(), BPPARAM, optim=FALSE)
    
    pdf('high_lambda_values_raw_counts.pdf')
    sapply(which(lambda(ods) > 1), function(i){
        plotExpressionRank(ods, normalized = FALSE, basePlot=TRUE, geneID=i,
                main=paste(rownames(ods)[i], 'lambda:', lambda(ods)[i]))
    })
    dev.off()
    optim
    optim <- TRUE
    res1 <- res
    res2 <- res
    i <- 402
    plotExpressionRank(ods, i, basePlot=TRUE)
    boxplot(lapply(1:ncol(res), function(x) res[,x]), log='y')
    plot(apply(res, 2, mean, trim=0.1))
    plot(colMedians(res))
    which(signif(colMeans(res), 6) == min(signif(colMeans(res), 6)))
}
