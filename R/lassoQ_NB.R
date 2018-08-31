
##
#https://de.mathworks.com/matlabcentral/fileexchange/34122-interface-to-a-lbfgs-solver?focused=5205992&tab=function

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
    fitls <- bplapply(seq_along(ods), lassoGeneCV2, D=D, b=b, k=k, sf=sf, H=H, theta=theta, 
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


absGradientsatZero <- function(fold, H, k, sf, theta, thetaC){
    k <- k[-fold]
    H <- H[-fold,]
    sf <- sf[-fold]
    b <- log(mean(k))
    abs(gradientD(c(b, rep(0, ncol(H))), H, k, sf, 1, theta, thetaC)[-1])
}

lassoGeneCV2 <- function(i, D, b, k, theta, lambda, H , sf, folds, optim, setPar=TRUE, debug=FALSE, ...){
    #pari3 <- c(b[i], D[i,])
    ki <- k[,i]
    nlambda <- 20
    
    # matrix to store results.
    res <- matrix(0, nrow=length(folds), ncol=nlambda)
    
    # estimate initial theta and mu parameters. 
    thetai <- estTheta(1, t(matrix(ki)), t(sf*mean(ki/sf)), H, 1, c(0.1, 200), NULL, negLogLikelihoodTheta)
    thetai <- thetai$minimum
    #initpar <- c(log(mean(ki)), numeric(length=ncol(H)))
    initpar <- c(b[i], D[i,])
    
    absGrad <- sapply(folds, absGradientsatZero, H=H, k=k, sf=sf, theta=thetai, thetaC=1)
    maxLambda <- max(absGrad)
    lambda <- c(exp(seq(log(maxLambda), -30, length.out=19)), 0)
    
    for(f in seq_along(folds)){
        
        pari <- initpar
        thetal <- thetai
        
        for(l in seq_along(lambda)){
            
            testSet <- folds[[f]]
            lam <- lambda[l]
            
            pari <- initpar
            #pari <- pari + rnorm(length(pari),0, 1E-8)
            fitpar <- lbfgs(truncLogLiklihoodD, gradientD, pari,  k=ki[-testSet], H=H[-testSet,], theta=thetal, 
                    exclusionMask=1, sf=sf[-testSet], thetaC=1, orthantwise_c=lam, epsilon = 1e-9,
                    orthantwise_start=1,orthantwise_end = length(pari), invisible=1, max_linesearch=50)
            pari <- fitpar$par
            
            mu <- predictGeneMu(pari, k[-testSet], H[-testSet,], sf[-testSet])
            thetal <- estTheta(1, t(matrix(ki[-testSet])), 
                    mu[,1], H[-testSet,], 1, 
                    c(0.1, 200), NULL, negLogLikelihoodTheta)$minimum
            
            fitpar <- lbfgs(truncLogLiklihoodD, gradientD, pari,  k=ki[-testSet], H=H[-testSet,], theta=thetal, 
                            exclusionMask=1, sf=sf[-testSet], thetaC=1, orthantwise_c=lam, epsilon = 1e-9,
                            orthantwise_start=2,orthantwise_end = length(pari), invisible=1, max_linesearch=50)
            pari <- fitpar$par
            
            # fitpar <- optim(pari, fn=truncLogLiklihoodDLasso, gr=gradientDLasso,
            #                 k=ki[-testSet], H=H[-testSet,], theta=thetal, exclusionMask=1, lambda=0,
            #                 sf=sf[-testSet], thetaC=1, lower=-100, upper=100, method='L-BFGS')

            mu <- predictGeneMu(pari, k[-testSet], H[-testSet,], sf[-testSet])
            thetal <- estTheta(1, t(matrix(ki[-testSet])), 
                   mu[,1], H[-testSet,], 1, 
                   c(0.1, 200), NULL, negLogLikelihoodTheta)$minimum
            
            
            #res[f, l] <- lossDtrunc(fitpar$par, ki[testSet], H[testSet,], sf[testSet], thetai, thetaC = 1)
            res[f, l] <- lossD(pari, ki[testSet], H[testSet,], sf[testSet], thetal)
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
    list('par'=c(fitpar$par), 'loss'=lossDtrunc(fitpar$par, ki[testSet], H[testSet,], sf[testSet], thetai, thetaC = 1))
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
