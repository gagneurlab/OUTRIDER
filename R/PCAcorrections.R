#'
#' the PCA implementation
#'
autoCorrectPCA <- function(ods, q, trim=0){
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- apply(x0, 2, mean, trim)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs=q) 
    pc  <- loadings(pca)
    
    # set the matrices
    E(ods) <- as.vector(pc)
    D(ods) <- as.vector(pc)
    b(ods) <- numeric(ncol(k))
    
    # add correction factor to object
    correctionFactors <- t(predictED(ods=ods))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    
    # add it to the object
    metadata(ods)[['dim']] <- dim(ods)
    
    validObject(ods)
    return(ods)
}

