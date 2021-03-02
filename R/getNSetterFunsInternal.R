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
    x0 <- transformed(ods)
    
    # compute per feature centered values
    b <- rowMeans(x0)
    x <- t(x0 - b)
    
    # add covariates as one-hot-encoded if requested
    if(!is.null(covariates(ods))){
        cov <- getCovariatesOneHotEncoded(ods)
        x <- cbind(x, cov)
    }
    
    return(x)
}

H <- function(ods){
    H <- x(ods) %*% E(ods)
    
    if(!is.null(covariates(ods))){
        H <- cbind(H, getCovariatesOneHotEncoded(ods))
    }
    
    return(H)
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
        assay(ods, 'replacedTrueCounts', withDimnames=FALSE) <- value
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

transformed <- function(ods){
    if(!("transformed" %in% assayNames(ods))){
        stop("Assay 'transformed' does not exist. Calculate it first with ", 
            "ods <- preprocess(ods)")
    }
    return(assay(ods, "transformed"))
}

`transformed<-` <- function(ods, value, ...){
    if(!is.matrix(value)){
        value <- matrix(value, nrow=nrow(ods))
    }
    assay(ods, "transformed", ...) <- value
    validObject(ods)
    return(ods)
}

getReverseTransformed <- function(x_trans, sf, revtransFUN){
    revtransFUN <- match.fun(revtransFUN)
    x <- revtransFUN(x_trans)
    x <- x * sf
    return(x)
}

covariates <- function(ods){
    return(metadata(ods)[['covariates']])
}

`covariates<-` <- function(ods, value){
    if(!all(value %in% colnames(colData(ods)))){
        stop("covariate has to be a column in colData(ods).")
    }
    metadata(ods)[['covariates']] <- value
    return(ods)
}

variability <- function(ods){
    if(profile(ods) == "outrider"){
        return(theta(ods))
    }
    else{
        return(rowSds(preprocessed(ods), na.rm=TRUE))
    }
}

