#' Coding conventions:
#' 
#' metadata(ods)[['E']]
#' metadata(ods)[['D']]
#' metadata(ods)[['b']]
#'
testing <- function(){
    devtools::load_all("~/projects/OUTRIDER")
    source("./../scared-analysis/src/r/helperFunction/simulateCountData.R")
    library(BBmisc)
    library(pcaMethods)
    set.seed(123)
    q <- 5
    simData <- createSimulationData(encDim=q, theta=c(mu=30, size=2), lognorm=FALSE)
    ods <- OutriderDataSet(countData=simData$k)
    assay(ods, 'inj_mask') <- simData$index
    assay(ods, 'trueCounts') <- round(simData$mu)
    mcols(ods)[['trueTheta']] <- simData$theta
    ods <- estimateSizeFactors(ods)
    theta <- 25
    control <- list(maxit=30)
    BPPARAM=MulticoreParam(20, progressbar=TRUE)
} 

autoCorrectED <- function(ods, q, theta=25, control=list(), BPPARAM=bpparam(), ...){
    
    # Check input
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    if(!'maxit' %in% names(control)){
        control$maxit <- 10
    }
    if(isScalarValue(theta)){
        theta <- rep(theta, nrow(ods))
    }
    
    k <- t(counts(ods, normalized=FALSE))
    sf <- sizeFactors(ods)
    
    # initialize W using PCA and bias as zeros.
    ods <- initED(ods, q=q, theta=theta)
    
    # check initial loss
    print(paste0('Initial PCA loss: ',
            lossED(getw(ods), k, getx(ods), sf, getb(ods), theta)))
    
    # optimize log likelihood
    t1 <- Sys.time()
    i <- 1
    for(i in 1:10){
        t2 <- Sys.time()
        
        theta <- updateTheta(ods, theta, BPPARAM)
        ods <- updateD(ods, theta, control, BPPARAM)
        ods <- updateE(ods, theta, control, BPPARAM)
        
        print(paste0('nb-PCA loss: ',
                lossED(getw(ods), k, getx(ods), s, getb(ods), theta)))
        print(Sys.time() - t2)
    }
    
    print(Sys.time() - t1)
    print(
        paste0('nb-PCA loss: ',
               lossED(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictED(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}

initED <- function(ods, q, theta, usePCA=TRUE){
    pc <- NULL
    if(isTRUE(usePCA)){
        pca <- pca(getx(ods), nPcs=q)
        pc  <- loadings(pca)
    } else {
        fit <- autoCorrectFit(getw(ods), loss2, lossGrad2, k, getx(ods), s,
                xbar, theta, control=list(factr <- 1E9, maxit=20))
        pc <- fit$par[1:prod(dim(ods))]
    }
    
    ods <- setD(ods, pc)
    ods <- setE(ods, pc)
    ods <- setb(ods, rowMeans(log(counts(ods) + 1)))
    
    return(ods)
}

updateTheta <- function(ods, theta, BPPARAM=bpparam()){
    return(theta)    
}

lossED <- function(w, k, x, sf, xbar, theta, minMu=0.01, ...){
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    ## encoding 
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    E <- W[,1:(ncol(W)/2)]
    D <- W[,(ncol(W)/2+1):ncol(W)]
    #theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    y <- t(t(x%*%E %*% t(D)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- minMu + sf*exp(y)
    
    ## log likelihood 
    ll <- dnbinom(t(k), mu=t(y_exp), size=theta, log=TRUE)
    - mean( ll )
}


lossGradED <- function(w, k, x, s, xbar, theta, ...){
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    E <- W[,1:(ncol(W)/2)]
    D <- W[,(ncol(W)/2+1):ncol(W)]
    theta <- matrix(theta, ncol=ncol(k), nrow=nrow(k), byrow=TRUE)
    
    # dW:
    t1 <- t(x) %*% (k %*% D)
    #t1 <- armaMatMultAtBC(x, k, W)
    t2 <- t(k) %*% (x %*% E)
    #t2 <- armaMatMultAtBC(k, x, W)
    y <- t(t(x%*%E %*% t(D)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)
    kt <- (k + theta)*y_exp/(y_exp+theta)
    
    t3 <- t(x) %*% (kt %*% D)
    t4 <- t(kt) %*% (x %*% E)
    
    dE <- (-t1  + t3)/prod(dim(k))
    dD <- (- t2 + t4)/prod(dim(k))
    #db:
    db <- colSums(kt-k)/prod(dim(k))
    
    if(!all(is.finite(dE)) | !all(is.finite(dD))){
        browser()
    }
    return(c(dE,dD, db))
}

predictED <- function(E, D, b, x, sf, ods, minMu=0.01){
    if(!missing(ods)){
        E <- getE(ods)
        D <- getD(ods)
        b <- getb(ods)
        x <- getx(ods)
        s <- sizeFactors(ods)
    }
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * (minMu + exp(y))
    y_exp
}
