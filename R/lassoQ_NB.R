


### Lasso CV to identify lambda.




updateLambda <- function(ods, control, BPPARAM){
    D <- D(ods)
    b <- b(ods)
    H <- H(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    mask <- exclusionMask(ods)
    theta <- theta(ods)
    thetaC <- thetaCorrection(ods)
    
    lambda = exp(seq(-5, 5, length.out = 20))
    
    fitls <- bplapply(1:nrow(ods), geneCV, D=D, b=b, k=k, sf=sf, H=H, 
                      theta=theta, mask=mask, control=control, thetaC=thetaC, 
                      lambda=lambda, BPPARAM=BPPARAM)
    
    # update D and bias terms
    # parMat <- sapply(fitls, '[[', 'par')
    # mcols(ods)['FitDMessage'] <- sapply(fitls, '[[', 'message')
    # print(table(mcols(ods)[,'FitDMessage']))
    # mcols(ods)[,'NumConvergedD'] <- mcols(ods)[,'NumConvergedD'] + grepl(
    #     "CONVERGENCE: REL_REDUCTION_OF_F .. FACTR.EPSMCH", 
    #     mcols(ods)[,'FitDMessage'])
    # b(ods) <- parMat[1,]
    # D(ods) <- t(parMat)[,-1]
    # 
    # metadata(ods)[['Dfits']] <- fitls
    lambda(ods) <- sapply(fitls, '[[', 'lambda')
    return(ods)
}


geneCV <- function(i, D, b, k, theta, mask, lambda, H , sf, ...){
    pari <- c(b[i], D[i,])
    ki <- k[,i]
    thetai <- theta[i]
    
    
    folds <- cut(seq(length(ki)),breaks=10,labels=FALSE)
    fold = 1:10
    ## lambda 
    res <- matrix(0, length(fold), length(lambda))
    for(l in seq(lambda)){
        for(f in fold){
            res[f,l] <- CVrun(lambda[l],f, folds=folds, ki=ki, H=H, theta=thetai, pari=pari, sf=sf)
        }
    }
    cmean = colMeans(res)
    which.min(cmean)
    return(c('means'=colMeans(res), 'min'= which.min(cmean), 'lambda'=lambda[which.min(cmean)]))
    
}


CVrun <- function(lambda, fold, folds, ki, H, theta, pari, sf){
    
    
    fitpar <- optim(pari, fn=truncLogLiklihoodDLasso, gr=gradientDLasso, 
                 k=ki[fold!=folds], H=H[fold!=folds,], theta=theta, exclusionMask=1, lambda=lambda,
                 sf=sf[fold!=folds], thetaC=1, lower=-100, upper=100, method='L-BFGS')
    lossD(fitpar$par, ki[fold==folds], H[fold==folds,], sf[fold==folds], theta)

}

lossD <- function(par, k, H, sf, theta){
    b <- par[1]
    d <- par[-1]
    
    y <- H %*% d + b
    yexp <- sf * exp(y)
    
    ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    
    # if(!is.finite(ll) & debugMyCode){
    #     browser()
    # }
    
    return(-ll)
}



