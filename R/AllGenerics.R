#' @rdname results
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))


counts.replace.OutriderDataSet <- function(object, value){
    assays(object)[["counts"]] <- value
    validObject(object)
    object
}

counts.OutriderDataSet <- function(object, normalized=FALSE, offset = 0) {
    
    cnts <- assays(object)[["counts"]]
    if(!normalized) {
        return(cnts)
    }
    if(!is.null(normalizationFactors(object))) {
        return((cnts + offset)/ normalizationFactors(object)*
                   rowMeans(normalizationFactors(object)))
    }
    if(is.null(sizeFactors(object)) || any(is.na(sizeFactors(object)))) {
        stop(paste("first calculate size factors, add normalizationFactors,",
                "or set normalized=FALSE"))
    }
    return( t( t( cnts + offset) / sizeFactors(object) * 
                   mean(sizeFactors(object)) ) )
}

setMethod("counts", signature(object="OutriderDataSet"),
          counts.OutriderDataSet)




#' 
#' Accessors for the 'counts' slot of a OutriderDataSet object.
#' 
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (gene or the like), and one
#' column for each sample. 
#'
#' By default this function returns the raw counts.
#' If conrols/normalizations are computed the normalized counts can be returned 
#' using normalized = TRUE.
#' The offset parameter can be used to add a pseudocount to the count before 
#' dividing by the normalization. This can be usefull when the log(counts) 
#' should be computed and in case the controll values are in the same oder of 
#' magnited as the counts.
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,OutriderDataSet-method 
#'         counts<-,OutriderDataSet,matrix-method
#'
#' @param object OutriderDataSet
#' @param normalized TRUE/FALSE whether counts should be normalized
#' @param offset pseudocount offset by default 0.
#'
#' @seealso \code{\link{sizeFactors}}, \code{\link{normalizationFactors}}
#'
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet(n=111, m=10)
#' counts(ods)[1:10,1:10]
#'
#' ods <- estimateSizeFactors(ods)
#' counts(ods, normalized=TRUE)[1:10,1:10]
#'
#' @export counts
setMethod("counts", signature(object="OutriderDataSet"), counts.OutriderDataSet)

#' @rdname counts
#' @export "counts<-"
setReplaceMethod("counts", signature(object="OutriderDataSet", value="matrix"),
        counts.replace.OutriderDataSet)   


normalizationFactors.OutriderDataSet <- function(object) {
    if (!"normalizationFactors" %in% assayNames(object)){
        return(NULL)
    }
    assays(object)[["normalizationFactors"]]
}

normFactors.replace.OutriderDataSet <- function(object, value, replace=TRUE) {
    # enforce same dimnames and matrix type
    if(!is.matrix(value)){
        value <- as.matrix(value)
    }
    dimnames(value) <- dimnames(object)
    
    # sanity checks
    stopifnot(!any(is.na(value)))
    stopifnot(all(is.finite(value)))
    stopifnot(all(value > 0))
    
    
    # multiply the new values with existing ones if not otherwise requested
    normF <- normalizationFactors(object)
    sizeF <- sizeFactors(object)
    if(replace == FALSE & !(is.null(normF) & is.null(sizeF))){
        if(!is.null(normF)){
            value <- value * normF
        } else {
            value <- t(t(value) * sizeF)
        }
    }
    
    # set the values and check the object
    assays(object)[["normalizationFactors"]] <- value
    validObject(object)
    object
}

#' 
#' Accessor functions for the normalization factors in a OutriderDataSet
#' object.
#'
#' @seealso DESeq2::normalizationFactors
#' @docType methods
#' @name normalizationFactors
#' @rdname normalizationFactors
#' @aliases normalizationFactors normalizationFactors,OutriderDataSet-method 
#'         normalizationFactors<-,OutriderDataSet,matrix-method 
#'         normalizationFactors<-,OutriderDataSet,DataFrame-method 
#'         normalizationFactors<-,OutriderDataSet,NULL-method
#' @param object a \code{OutriderDataSet} object.
#' @param value the matrix of normalization factors
#' @param replace if old values are present values are replcaed. If set to false
#'                old and new values are multiplied.      
#' @param ... additional arguments
#' @examples
#'
#' ods <- makeExampleOutriderDataSet(n=111, m=10)
#'
#' normFactors <- matrix(runif(nrow(ods)*ncol(ods),0.5,1.5),
#'     ncol=ncol(ods),nrow=nrow(ods))
#'
#' # the normalization factors matrix should not have 0's in it
#' # it should have geometric mean near 1 for each row
#' normFactorsRM <- normFactors / exp(rowMeans(log(normFactors)))
#' normalizationFactors(ods) <- normFactorsRM
#' normalizationFactors(ods)[1:10,1:10]
#' 
#' normalizationFactors(ods) <- NULL
#' ods <- estimateSizeFactors(ods)
#' normalizationFactors(ods, replace=FALSE) <- normFactors
#' all(normalizationFactors(ods) == t(sizeFactors(ods) * t(normFactors)))
#' 
#' @export normalizationFactors
setMethod("normalizationFactors", signature(object="OutriderDataSet"),
        normalizationFactors.OutriderDataSet)

#' @rdname normalizationFactors
#' @export "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="OutriderDataSet", 
        value="matrix"), normFactors.replace.OutriderDataSet)

#' @rdname normalizationFactors
#' @export "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="OutriderDataSet", 
        value="DataFrame"), normFactors.replace.OutriderDataSet)

#' @rdname normalizationFactors
#' @export "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="OutriderDataSet", 
        value="data.frame"), normFactors.replace.OutriderDataSet)

#' @rdname normalizationFactors
#' @export "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="OutriderDataSet", 
        value="NULL"), function(object, value) {
                assays(object)[["normalizationFactors"]] <- NULL
                return(object)})

