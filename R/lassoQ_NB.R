#'
#' Lasso CV to identify lambda.
#' 
#' @noRd
updateLambda <- function(ods, nFolds=10, control, BPPARAM){
    D <- D(ods)
    b <- b(ods)
    H <- H(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    theta <- theta(ods)
    
    lambda = exp(c(-10, -7, seq(-5, 5, length.out=18)))
    t1 <- Sys.time()
    fitls <- bplapply(seq_along(ods), lassoGeneCV, D=D, b=b, k=k, sf=sf, H=H, theta=theta, 
            lambda=lambda, nFolds=nFolds, control=control, BPPARAM=BPPARAM)
    t2 <- Sys.time()
    
    lambda(ods) <- sapply(fitls, '[[', 'lambda')
    
    print(paste('Lasso lambda fit time: ', t2 - t1))
    print(summary(lambda(ods)))
    
    return(ods)
}


lassoGeneCV <- function(i, D, b, k, theta, lambda, H , sf, nFolds=20, ...){
    pari <- c(b[i], D[i,])
    ki <- k[,i]
    thetai <- theta[i]
    
    folds <- BBmisc::chunk(sample(seq_along(ki)), n.chunks=nFolds)
    
    res <- sapply(lambda, function(l){
        sapply(seq_len(nFolds), function(f){
            lassoSingleCVrun(l, testSet=folds[[f]], ki=ki, H=H, theta=thetai, 
                    pari=pari, sf=sf)
        })
    })
    
    cmean = colMeans(res)
    ans <- c('means'=colMeans(res), 'minIdx'= which.min(cmean), 
      'lambda'=lambda[which.min(cmean)])
    return(ans)
}

lassoSingleCVrun <- function(l, testSet, ki, H, thetai, pari, sf){
    fitpar <- optim(pari, fn=truncLogLiklihoodDLasso, gr=gradientDLasso,
            k=ki[-testSet], H=H[-testSet,], theta=thetai, exclusionMask=1, 
            lambda=l, sf=sf[-testSet], thetaC=1, lower=-100, upper=100, 
            method='L-BFGS')
    lossD(fitpar$par, ki[testSet], H[testSet,], sf[testSet], thetai)
}

if(FALSE){
    pdf('high_lambda_values_raw_counts.pdf')
    sapply(which(lambda(ods) > 1), function(i){
        plotExpressionRank(ods, normalized = FALSE, basePlot=TRUE, geneID=i,
                main=paste(rownames(ods)[i], 'lambda:', lambda(ods)[i]))
    })
    dev.off()
}