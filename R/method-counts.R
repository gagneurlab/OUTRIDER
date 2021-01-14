
counts.replace.OutriderDataSet <- function(object, ..., value){
    mode(value) <- "integer"
    assay(object, "counts", ...) <- value

    validObject(object)
    object
}

counts.OutriderDataSet <- function(object, normalized=FALSE, minE=0.5, ...){
    cnts <- assay(object, "counts", ...)
    
    # raw counts
    if(!normalized) {
        return(cnts)
    }
    
    # normalized by normalization factors
    if(!is.null(normalizationFactors(object))) {
        if(!"expectedLogGeomMean" %in% colnames(mcols(object))){
            stop("The expectedLogGeomMean is missing in mcols(ods). ",
                    "Did you set the normaliationFactors by hand? ",
                    "Please use normalizationFactors(object) <- values.")
        }
        
        # use cached expected log geom mean values
        eMat <- pmax(normalizationFactors(object), minE)
        eLGM <- mcols(object)[["expectedLogGeomMean"]]
        return(cnts/eMat * eLGM)
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
#' @param object An \code{\link{OutriderDataSet}} object
#' @param normalized TRUE/FALSE whether counts should be normalized
#' @param value An integer matrix containing the counts
#' @param minE The minimal expected count, defaults to 0.5, to be used in 
#'          computing the expected log geom mean.
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


normalizationFactors.Outrider2DataSet <- function(object, ...) {
    if (!"normalizationFactors" %in% assayNames(object)){
        return(NULL)
    }
    assay(object, "normalizationFactors", ...)
}

normFactors.replace.OutriderDataSet <- function(object, minE=0.5, ..., value) {
    # enforce same dimnames and matrix type
    if(!is.matrix(value)){
        value <- as.matrix(value)
    }
    dimnames(value) <- dimnames(object)
    
    # sanity checks
    if(modelParams(object)$distribution == "negative binomial"){
        stopifnot(!any(is.na(value)))
        stopifnot(all(is.finite(value)))
        stopifnot(all(value > 0))
    }
    
    # set the values and check the object
    assay(object, "normalizationFactors", ..., withDimnames=FALSE) <- value
    
    # compute the expected log geom mean values so we can cache them
    mcols(object)[["expectedLogGeomMean"]] <- exp(
        rowMeans2(log(pmax(value, minE)), na.rm=TRUE))
    
    # validate and return object
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
#' @docType methods
#' @name normalizationFactors
#' @rdname normalizationFactors
#' @aliases normalizationFactors normalizationFactors,OutriderDataSet-method 
#'         normalizationFactors<-,OutriderDataSet,matrix-method 
#'         normalizationFactors<-,OutriderDataSet,DataFrame-method 
#'         normalizationFactors<-,OutriderDataSet,NULL-method
#' 
#' @inheritParams counts
#' @param value The matrix of normalization factors
#' @return A numeric matrix containing the normalization factors or the 
#'             OutriderDataSet object with an updated 
#'             \code{normalizationFactors} assay.
#' 
#' @seealso \code{\link{sizeFactors}} 
#'          \code{\link[DESeq2]{normalizationFactors}}
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
setMethod("normalizationFactors", signature(object="Outrider2DataSet"),
        normalizationFactors.Outrider2DataSet)

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
setReplaceMethod("normalizationFactors", signature(object="Outrider2DataSet", 
        value="matrix"), normFactors.replace.OutriderDataSet)

#' @rdname normalizationFactors
#' @export "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="Outrider2DataSet", 
        value="DataFrame"), normFactors.replace.OutriderDataSet)

#' @rdname normalizationFactors
#' @export "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="Outrider2DataSet", 
        value="data.frame"), normFactors.replace.OutriderDataSet)

#' @rdname normalizationFactors
#' @export "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="Outrider2DataSet", 
        value="NULL"), function(object, value) {
                assay(object, "normalizationFactors") <- NULL
                return(object)})

