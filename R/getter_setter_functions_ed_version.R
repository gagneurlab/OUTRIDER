
getx <- function(ods){
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    b <- colMeans(x0)
    x <- t(t(x0) - b)
    
    return(x)
}

setD <- function(ods, D){
    if(!is.matrix(D)){
        D <- matrix(D, nrow=nrow(ods))
    }
    metadata(ods)[['D']] <- D
    return(ods)
}

getD <- function(ods){
    return(metadata(ods)[['D']])
}

setE <- function(ods, E){
    if(!is.matrix(E)){
        E <- matrix(E, nrow=nrow(ods))
    }
    metadata(ods)[['E']] <- E
    return(ods)
}

getE <- function(ods){
    return(metadata(ods)[['E']])
}

setb <- function(ods, b){
    mcols(ods)[['b']] <- b
    return(ods)
}

getb <- function(ods){
    return(mcols(ods)[['b']])
}

getw <- function(ods){
    return(c(as.vector(getE(ods)), as.vector(getD(ods)), getb(ods)))
}

getTrueCounts <- function(ods){
    if('replacedTrueCounts' %in% assayNames(ods)){
        return(assay(ods, 'replacedTrueCounts'))
    }
    return(counts(ods))
}

setTrueCounts <- function(ods, cts){
    if(!'replacedTrueCounts' %in% assayNames(ods)){
        assay(ods, 'replacedTrueCounts') <- cts
    }
    return(ods)
}

getExclusionMask <- function(ods){
    if('exclusionMask' %in% assayNames(ods)){
        return(assay(ods, 'exclusionMask'))
    }
    return(matrix(1, ncol=ncol(ods), nrow=nrow(ods)))
}

setExclusionMask <- function(ods, mask){
    assay(ods, 'exclusionMask') <- mask
    return(ods)
}