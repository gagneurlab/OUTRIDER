autoCorrectPCA <- function(ods, q){
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w <- c(as.vector(pc), numeric(ncol(k)))
    
    correctionFactors <- t(predictC(w, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}

autoCorrectPCAtrimmedMean <- function(ods, q){
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    # use 5% trimmed mean to center x. 
    xbar <- apply(x0, 2, mean, 0.05)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w <- c(as.vector(pc), numeric(ncol(k)))
    
    correctionFactors <- t(predictC(w, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}


autoCorrectPCACooks <- function(ods, q, predictOnReplacedCounts=FALSE){
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    rep_k <- replaceOutliersCooks(k, BPPARAM=bpparam())
    k_no <- rep_k$cts
    
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w <- c(as.vector(pc), numeric(ncol(k)))
    
    if(isTRUE(predictOnReplacedCounts)){
        correctionFactors <- t(predictC(w, k_no, s, xbar))
    }else{
        correctionFactors <- t(predictC(w, k, s, xbar))
    }
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}

# trimming the mean and DESeq Cooks replacemet gives bad results, especially
# low expression outliers.
autoCorrectPCAtrimmedMeanCooks <- function(ods, q, predictOnReplacedCounts=FALSE){
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    rep_k <- replaceOutliersCooks(k, BPPARAM=bpparam())
    k_no <- rep_k$cts
    
    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
    # use 5% trimmed mean to center x. 
    xbar <- apply(x0, 2, mean, 0.05)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w <- c(as.vector(pc), numeric(ncol(k)))
    
    if(isTRUE(predictOnReplacedCounts)){
        correctionFactors <- t(predictC(w, k_no, s, xbar))
    }else{
        correctionFactors <- t(predictC(w, k, s, xbar))
    }
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}


autoCorrectPCAstandardized <- function(ods, q){
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    xSd <- colSds(x0)
    x <- t((t(x0) - xbar)/xSd)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w <- c(as.vector(pc), numeric(ncol(k)))
    
    
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    y <- t(t(x%*%W %*% t(W)) * xSd +b + xbar)
    correctionFactors <- t(s*exp(y))
    #correctionFactors <- t(predictC(w, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}
