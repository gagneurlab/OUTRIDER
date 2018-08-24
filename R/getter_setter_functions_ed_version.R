x <- function(ods){
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    b <- colMeans(x0)
    x <- t(t(x0) - b)
    
    return(x)
}

H <- function(ods){
    x(ods) %*% E(ods)
}

`D<-` <- function(ods, value){
    if(!is.matrix(value)){
        value <- matrix(value, nrow=nrow(ods))
    }
    metadata(ods)[['D']] <- value
    return(ods)
}

D <- function(ods){
    return(metadata(ods)[['D']])
}

`E<-` <- function(ods, value){
    if(!is.matrix(value)){
        value <- matrix(value, nrow=nrow(ods))
    }
    metadata(ods)[['E']] <- value
    return(ods)
}

E <- function(ods){
    return(metadata(ods)[['E']])
}

`b<-` <- function(ods, value){
    mcols(ods)[['b']] <- value
    return(ods)
}

b <- function(ods){
    return(mcols(ods)[['b']])
}

predictC <- function(ods){
    predictMatC(x(ods), E(ods), D(ods), b(ods), sizeFactors(ods))
}

predictY <- function(ods){
    predictMatY(x(ods), E(ods), D(ods), b(ods))
}

getw <- function(ods){
    return(c(as.vector(getE(ods)), as.vector(getD(ods)), getb(ods)))
}

trueCounts <- function(ods){
    if('replacedTrueCounts' %in% assayNames(ods)){
        return(assay(ods, 'replacedTrueCounts'))
    }
    return(counts(ods))
}

`trueCounts<-` <- function(ods, value){
    if(!'replacedTrueCounts' %in% assayNames(ods)){
        assay(ods, 'replacedTrueCounts') <- value
    }
    return(ods)
}

thetaCorrection <- function(ods){
    if(!"thetaCorrection" %in% colnames(colData(ods))){
        warning('thetaFactors are not computed. If this intended you can ', 
                'ignore this message by setting them to 1. Otherwise please ',
                'fit the autoencoder first.')
        return(rep(1, ncol(ods)))
    }
    return(colData(ods)[,'thetaCorrection'])
}

`thetaCorrection<-` <- function(ods, value){
    if(isScalarNumeric(value)){
        value <- rep(value, ncol(ods))
    }
    colData(ods)[,'thetaCorrection'] <- value
    return(ods)
}

lambda <- function(ods){
    if(!"lambda" %in% colnames(mcols(ods))){
        warning('TODO Replace by apropriate warning.')
        return(rep(0, nrow(ods)))
    }
    return(mcols(ods)[,'lambda'])
}

`lambda<-` <- function(ods, value){
    if(isScalarNumeric(value)){
        value <- rep(value, nrow(ods))
    }
    mcols(ods)[,'lambda'] <- value
    return(ods)
}


exclusionMask <- function(ods){
    if('exclusionMask' %in% assayNames(ods)){
        return(assay(ods, 'exclusionMask'))
    }
    return(matrix(1, ncol=ncol(ods), nrow=nrow(ods)))
}

`exclusionMask<-` <- function(ods, value){
    if(isScalarNumeric(value)){
        value <- matrix(value, ncol=ncol(ods), nrow=nrow(ods))
    }
    assay(ods, 'exclusionMask') <- value
    return(ods)
}