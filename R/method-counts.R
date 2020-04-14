
counts.replace.OutriderDataSet <- function(object, ..., value){
    mode(value) <- "integer"
    assay(object, "counts", ...) <- value

    validObject(object)
    object
}

counts.OutriderDataSet <- function(object, normalized=FALSE, minE=0.5){
    cnts <- assay(object, "counts")
    
    # raw counts
    if(!normalized) {
        return(cnts)
    }
    
    # normalized by normalization factors
    if(!is.null(normalizationFactors(object))) {
        E <- t(apply(normalizationFactors(object), 1, pmax, minE))
        return(cnts/E * exp(rowMeans(log(E))))
    }
    
    # normalization by sizeFactors
    if(is.null(sizeFactors(object)) || any(is.na(sizeFactors(object)))) {
        stop(paste("first calculate size factors, add normalizationFactors,",
                "or set normalized=FALSE"))
    }
    return(t(t(cnts) / sizeFactors(object)))
}

#' 
#' Accessors for the 'counts' slot of an OutriderDataSet object.
#' 
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (eg.: gene), and one
#' column for each sample. 
#'
#' By default this function returns the raw counts. If conrol factors are
#' computed or provided the normalized counts can be returned using 
#' normalized = TRUE. The offset parameter can be used to add a pseudocount
#' to the count before dividing by the normalization. This can be usefull 
#' when the log(counts) are computed and in case the controll values are in 
#' the same oder of magnited as the counts.
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,OutriderDataSet-method 
#'         counts<-,OutriderDataSet,matrix-method
#'
#' @param object OutriderDataSet
#' @param normalized TRUE/FALSE whether counts should be normalized
#' @param value An integer matrix containing the counts
#' @param minE minimal expected count.
#' @param ... Further arguments are passed on to the underlying assay function
#' @return A matrix containing the counts
#' 
#' @seealso \code{\link{sizeFactors}}, \code{\link{normalizationFactors}}
#'
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet()
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
    assay(object, "normalizationFactors")
}

normFactors.replace.OutriderDataSet <- function(object, value) {
    # enforce same dimnames and matrix type
    if(!is.matrix(value)){
        value <- as.matrix(value)
    }
    dimnames(value) <- dimnames(object)
    
    # sanity checks
    stopifnot(!any(is.na(value)))
    stopifnot(all(is.finite(value)))
    stopifnot(all(value > 0))
    
    # set the values and check the object
    assay(object, "normalizationFactors", withDimnames=FALSE) <- value
    validObject(object)
    object
}

#' 
#' Accessor functions for the normalization factors in an OutriderDataSet
#' object.
#' 
#' To normalize raw count data normalization factors can be provided as
#' a matrix. When running \code{\link{controlForConfounders}} the normalization 
#' factors are stored within the OutriderDataset object. This normalization 
#' factors are then used to compute the normalized counts.
#'
#' @seealso DESeq2::normalizationFactors
#' @docType methods
#' @name normalizationFactors
#' @rdname normalizationFactors
#' @aliases normalizationFactors normalizationFactors,OutriderDataSet-method 
#'         normalizationFactors<-,OutriderDataSet,matrix-method 
#'         normalizationFactors<-,OutriderDataSet,DataFrame-method 
#'         normalizationFactors<-,OutriderDataSet,NULL-method
#'         
#' @param object An \code{OutriderDataSet} object.
#' @param value The matrix of normalization factors
#' @return A numeric matrix containing the normalization factors or the 
#'             OutriderDataSet object with an updated 
#'             \code{normalizationFactors} assay.
#' 
#' @seealso \code{\link{sizeFactors}}
#' 
#' @examples
#'
#' ods <- makeExampleOutriderDataSet()
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
#' normalizationFactors(ods) <- normFactors
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
                assay(object, "normalizationFactors") <- NULL
                return(object)})

