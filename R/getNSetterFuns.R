#' Getter/Setter functions
#' 
#' This is a collection of small accessor/setter functions for easy access to
#' the values within the OUTRIDER model.
#' 
#' @param ods,object An OutriderDataSet object.
#' @param ... Further arguments passed on to the underlying assay function.
#' @return A matrix or vector dependent on the type of data retrieved.
#' 
#' @name getter_setter_functions
#' @rdname getter_setter_functions
#' @aliases zScore, pValue, padj, theta, dispersion, getBestQ
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(10, 10)
#' ods <- OUTRIDER(ods)
#' 
#' zScore(ods)
#' pValue(ods)
#' padj(ods)
#' theta(ods)
#' theta(ods) == 1/dispersions(ods)
#' getBestQ(ods)
#' 
NULL

#' @rdname getter_setter_functions
#' @export
zScore <- function(ods){
    if(!'zScore' %in% assayNames(ods)){
        stop('Please compute first the Z-scores before retrieving them.')
    }
    assay(ods, 'zScore')
}

`zScore<-` <- function(ods, ..., value){
    stopifnot(is.matrix(value))
    stopifnot(dim(ods) == dim(value))
    assay(ods, 'zScore', ...) <- value
    return(ods)
}

#' @rdname getter_setter_functions
#' @export
pValue <- function(ods){
    if(!'pValue' %in% assayNames(ods)){
        stop('Please compute first the P-values before retrieving them.')
    }
    assay(ods, 'pValue')
}

`pValue<-` <- function(ods, ..., value){
    stopifnot(is.matrix(value))
    stopifnot(dim(ods) == dim(value))
    assay(ods, 'pValue', ...) <- value
    return(ods)
}

#' @rdname getter_setter_functions
#' @export padj
padj <- function(ods){
    if(!'padjust' %in% assayNames(ods)){
        stop('Please compute first the P-values before retrieving', '
                the adjusted ones.')
    }
    assay(ods, 'padjust')
}

`padj<-` <- function(ods, ..., value){
    stopifnot(is.matrix(value))
    stopifnot(dim(ods) == dim(value))
    assay(ods, 'padjust', ...) <- value
    return(ods)
}

#' @rdname getter_setter_functions
#' @export dispersions
setMethod("dispersions", signature(object="OutriderDataSet"), 
    function(object, ...){ 1/theta(object) })

#' @rdname getter_setter_functions
#' @export theta
theta <- function(ods){
    if(!'theta' %in% colnames(mcols(ods))){
        stop('Please fit first the autoencoder before retrieving thetas.')
    }
    mcols(ods)[['theta']]
}

`theta<-` <- function(ods, value){
    mcols(ods)[['theta']] <- value
    return(ods)
}


#' @rdname getter_setter_functions
#' @export
setMethod("modelParams", "Outrider2DataSet", function(object, paramName) {
    # modelParams <- list(distribution = slot(object, "distribution"),
    #                     transformation = slot(object, "transformation"),
    #                     preprocessing = slot(object, "preprocessing"),
    #                     fitModel = slot(object, "fitModel"))
    if(!(paramName %in% c("distribution", "transformation", "preprocessing", 
                          "fitModel"))){
        stop("Argument paramName needs to be one of distribution, ", 
             "transformation, preprocessing or fitModel.")
    }
    modelParam <- slot(object, paramName)
    return(modelParam)
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("modelParams", "Outrider2DataSet", function(object, paramName, value) {
    slot(object, paramName) <- value
    validObject(object)
    return(object)
})

#' @rdname getter_setter_functions
#' @export
setMethod("observed", "Outrider2DataSet", function(object, normalized=FALSE, 
                                                    minE=0) {
    
    if(!("observed" %in% assayNames(object))){
        stop("Assay 'observed' does not exist.")
    }
    obs <- assay(object, "observed")
    
    if(isFALSE(normalized)){
        return(obs)
    }
    else{
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
            return(obs/eMat * eLGM)
        }
        
        # normalization by sizeFactors
        if(is.null(sizeFactors(object)) || any(is.na(sizeFactors(object)))) {
            stop(paste("first calculate size factors, add normalizationFactors,",
                       "or set normalized=FALSE"))
        }
        return(t(t(obs) / sizeFactors(object)))
    }
})

#' @rdname getter_setter_functions
#' @export
setMethod("observed", "OutriderDataSet", function(object, normalized=FALSE) {
    return(counts(object, normalized=normalized))
})

#' @rdname getter_setter_functions
#' @export
setMethod("observed", "ProtriderDataSet", function(object, normalized=FALSE) {
    return(assay(object, "intensities"))
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("observed", "Outrider2DataSet", function(object, value) {
    assay(object, "observed") <- value
    validObject(object)
    return(object)
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("observed", "OutriderDataSet", function(object, value) {
    counts(object, normalized=FALSE) <- value
    return(object)
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("observed", "ProtriderDataSet", function(object, value) {
    assay(object, "intensities") <- value
    validObject(object)
    return(object)
})


#' @rdname getter_setter_functions
#' @export
setMethod("preprocessed", "Outrider2DataSet", function(object) {
    if(!("preprocessed" %in% assayNames(object))){
        return(raw(ods))
    }
    return(assay(object, "preprocessed"))
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("preprocessed", "Outrider2DataSet", function(object, value) {
    assay(object, "preprocessed") <- value
    validObject(object)
    return(object)
})


