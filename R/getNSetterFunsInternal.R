#' Sample exclusion
#' 
#' To exclude a sample from the fit process, one can use this function to mask 
#' specific samples. This can be used if replicates are present.
#' 
#' @param ods An OutriderDataSet object
#' @param value A logical vector of the length of the samples. If \code{TRUE},
#'             the corresponding sample will be excluded from the autoencoder
#'             fit.
#' @param aeMatrix If \code{TRUE}, it returns a 0/1 matrix for the 
#'             internal autoencoder functions in the form of feature x sample
#' @return The exclusion vector/matrix.
#' 
#' @name sampleExclusionMask
#' @rdname sampleExclusionMask
#' @aliases sampleExclusionMask, `sampleExclusionMask<-`
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' sampleExclusionMask(ods) <- sample(c(FALSE, TRUE), ncol(ods), replace=TRUE)
#' 
#' sampleExclusionMask(ods)
#' 
#' @export sampleExclusionMask
sampleExclusionMask <- function(ods, aeMatrix=FALSE){
    if('exclude' %in% colnames(colData(ods))){
        ans <- colData(ods)[['exclude']]
    } else {
        ans <- rep(FALSE, ncol(ods))
    }
    names(ans) <- colnames(ods)
    
    if(isTRUE(aeMatrix)){
        ans <- as.integer(vapply(ans, isFALSE, FALSE))
        ans <- matrix(ans, ncol=ncol(ods), nrow=nrow(ods), byrow=TRUE)
        colnames(ans) <- colnames(ods)
        rownames(ans) <- rownames(ods)
    }
    
    return(ans)
}

#' @rdname sampleExclusionMask
#' @export "sampleExclusionMask<-"
`sampleExclusionMask<-` <- function(ods, value){
    if(isScalarLogical(value)){
        value <- rep(value, ncol(ods))
    }
    colData(ods)[['exclude']] <- value
    return(ods)
}

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
        #warning('thetaFactors are not computed. If this intended you can ', 
        #        'ignore this message by setting them to 1. Otherwise please ',
        #        'fit the autoencoder first.')
        return(rep(1, ncol(ods)))
    }
    return(colData(ods)[,'thetaCorrection'])
}

`thetaCorrection<-` <- function(ods, value){
    if(isScalarNumeric(value)){
        value <- rep(value, ncol(ods))
    }
    colData(ods)[['thetaCorrection']] <- value
    return(ods)
}
