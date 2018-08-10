#' Coding conventions:
#' 
#' metadata(ods)[['E']]
#' metadata(ods)[['D']]
#' metadata(ods)[['b']]
#'
testing <- function(){
    odsg <- readRDS("../scared-analysis/Output/data/GTEx_not_sun_exposed_OutriderDONE.RDS")
    odsg <- readRDS('/s/project/scared/paper/revision/run2107/data/Skin_Not_Sun_Exposed_Suprapubic_Outrider.RDS')
    ods <- odsg[1:5000, 1:100]
    q <- 10
    q <- 49
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
    control <- list(maxit=40)
    BPPARAM=SerialParam()
    BPPARAM=MulticoreParam(4, progressbar=TRUE, stop.on.error = FALSE)
} 

autoCorrectED <- function(ods, q, theta=25, control=list(), BPPARAM=bpparam(), 
        minMu=0.01, ...){
    lossList <- list()
    lossList2 <- list()
    printTheta <- 10
    # Check input
    if(!'factr' %in% names(control)){
        control$factr <- 1E7
    }
    if(!'maxit' %in% names(control)){
        control$maxit <- 100
    }
    if(isScalarValue(theta)){
        theta <- rep(theta, nrow(ods))
    }
    
    k <- t(counts(ods, normalized=FALSE))
    sf <- sizeFactors(ods)
    
    # initialize W using PCA and bias as zeros.
    ods <- initED(ods, q=q, theta=theta, usePCA=TRUE)
    #theta <- updateTheta(ods, theta, BPPARAM)
    
    # check initial loss
    print(paste0('Initial PCA loss: ', lossED(ods, theta)))
    lossList[1] <- lossED(ods, theta)
    # optimize log likelihood
    if(!bpisup(BPPARAM)){
        bpstart(BPPARAM)
    }
    t1 <- Sys.time()
    for(i in 1:10){
        t2 <- Sys.time()
        # check initial loss
        # print(paste0(i, 'Initial PCA loss: ', lossED(ods, theta)))
        #if(i==1){
        # for(a in 1:5){
        #         ods <- updateD(ods, theta, control, BPPARAM)
        #         theta <- updateTheta(ods, theta, BPPARAM)
        #         print(paste0(i, 'initGLM loss: ', lossED(ods, theta)))
        #     #}
        # }
        ods <- updateD(ods, theta, control, BPPARAM)
        
        # check initial loss
        print(paste0(i, 'update D loss: ', lossED(ods, theta)))
        lossList[1+i*3-2] <- lossED(ods, theta)
        
        theta <- updateTheta(ods, theta, BPPARAM)
        
        # check initial loss
        print(paste0(i, 'Theta loss: ', lossED(ods, theta)))
        lossList[1+i*3-1] <- lossED(ods, theta)
        
        
        ods <- updateE(ods, theta, control, BPPARAM)
        
        print(paste0(i, 'update E loss: ', lossED(ods, theta)))
        lossList[1+i*3] <- lossED(ods, theta)
        lossList2[i] <- lossED(ods, theta=10)
        print(Sys.time() - t2)
    }
    
    print(Sys.time() - t1)
    print(paste0(i, 'Final nb-PCA loss: ',
            lossED(ods, theta)))
    
    bpstop(BPPARAM)
    
    correctionFactors <- t(predictED(ods=ods))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- getw(ods)
    metadata(ods)[['dim']] <- dim(ods)
    metadata(ods)[['loss']] <- lossList
    metadata(ods)[['lossT10']] <- lossList2
    mcols(ods)[['disp']] <- theta
    validObject(ods)
    
    return(ods)
}

initED <- function(ods, q, theta, usePCA=TRUE){
    pca <- pca(getx(ods), nPcs=q)
    pc  <- loadings(pca)
    
    if(isFALSE(usePCA)){
        pc[1:length(pc)] <- rnorm(length(pc), sd=sd(pc))
    }
    
    ods <- setD(ods, pc)
    ods <- setE(ods, pc)
    ods <- setb(ods, rowMeans(log(counts(ods) + 1)))
    
    # if(isFALSE(usePCA)){
    #     ods <- updateE(ods, theta, control, BPPARAM)
    #     ods <- setD(ods, getE(ods))
    # }
    
    return(ods)
}

updateTheta <- function(ods, theta, BPPARAM=bpparam()){
    normalizationFactors(ods, replace=TRUE) <- t(predictED(ods=ods))
    theta <- dispersions(fit(ods, BPPARAM=BPPARAM))
    return(theta)    
}

lossED <- function(ods, theta, minMu=0.01, ...){
    
    ## encoding 
    b <- getb(ods)
    E <- getE(ods)
    D <- getD(ods)
    x <- getx(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    
    y <- t(t(x %*% E %*% t(D)) + b)
    #y_exp <- pmin(1e7, sf * (minMu + exp(y)))
    y_exp <- sf * (minMu + exp(y))
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
        sf <- sizeFactors(ods)
    }
    
    y <- t(t(x %*% E %*% t(D)) + b)
    y_exp <- sf * (minMu + exp(y))
    y_exp
}
