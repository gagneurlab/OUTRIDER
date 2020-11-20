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
#' @aliases dispersions dispersions,OutriderDataSet-method
#' @seealso \code{\link[DESeq2]{estimateDispersions}}
#' @export
setMethod("dispersions", signature(object="OutriderDataSet"),
    function(object, ...){ 1/theta(object) })

setMethod("dispersions", signature(object="Outrider2DataSet"),
          function(object, ...){ 
                if(modelParams(object, "distribution") != "negative binomial"){
                      stop("Dispersion is only defined for negative binomial, ", 
                           "not for ", modelParams(object, "distribution"))
                } else{
                    1/theta(object) 
                }
})


#' @rdname getter_setter_functions
#' @export theta
theta <- function(ods){
    if(modelParams(ods, "distribution") != "negative binomial"){
        stop("theta is only defined for negative- binomial distribution, ", 
             "not for ", modelParams(object, "distribution"))
    }
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
    if(missing(paramName)){
        modelParams <- list(distribution = slot(object, "distribution"),
                            preprocessing = slot(object, "preprocessing"),
                            transformation = slot(object, "transformation"),
                            sf_norm = slot(object, "sf_norm"))
        return(modelParams)
    
    } else if(!(paramName %in% c("distribution", "transformation", 
                                    "preprocessing", "sf_norm"))){
        stop("Argument paramName needs to be one of distribution, sf_norm, ", 
             "transformation or preprocessing.")
    }
    modelParam <- slot(object, paramName)
    return(modelParam)
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("modelParams", "Outrider2DataSet", 
        function(object, value, paramName) {
    params <- c("distribution", "preprocessing", "transformation", "sf_norm")
    if(!missing(paramName)){
        slot(object, paramName) <- value
    } else if(any(params %in% names(value))){
        for(param in params[params %in% names(value)]){
            slot(object, param) <- value[[param]] 
        }
    } else{
        stop("Either specify the parameter name that shall be replaced or ", 
             "provide a named list with the parameters to replace.")
    }
    validObject(object)
    return(object)
})

observed.OUTRIDER2 <- function(object, normalized=FALSE, minE=0.5, ...){
    
    if(modelParams(object, "distribution") != "negative binomial"){
        minE=-Inf
    }
    
    if(!("observed" %in% assayNames(object))){
        if(!("counts" %in% assayNames(object))){
            stop("Assay 'observed' does not exist.")
        } else{
            obs <- assay(object, "counts")
        }
    } else{
        obs <- assay(object, "observed")
    }
    
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
            
            # use preprocessed observations
            obs <- preprocessed(object)
            
            # normalize
            if(modelParams(object, "distribution") == "negative binomial"){
                return(obs/eMat * eLGM)
            } else{
                norm <- (obs - eMat) + eLGM
                if(modelParams(object, "transformation") == "log" ||
                   modelParams(object, "preprocessing") == "log"){
                    norm[norm < 0] <- 0
                }
                return(norm)
                
            }
            
            
        }
        
        # normalization by sizeFactors
        if(is.null(sizeFactors(object)) || any(is.na(sizeFactors(object)))) {
            stop(paste("first calculate size factors, add normalizationFactors,",
                       "or set normalized=FALSE"))
        }
        return(t(t(obs) / sizeFactors(object)))
    }
}

#' @rdname getter_setter_functions
#' @export
setMethod("observed", "Outrider2DataSet", observed.OUTRIDER2)

#' #' @rdname getter_setter_functions
#' #' @export
#' setMethod("observed", "OutriderDataSet",
#'     function(object, normalized=FALSE, minE=0.5, ...) {
#'         observed.OUTRIDER2(object, normalized=normalized, minE=minE, ...)
#' })

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("observed", "Outrider2DataSet", function(object, value, ...) {
    assay(object, "observed", ...) <- value
    validObject(object)
    return(object)
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("observed", "OutriderDataSet", function(object, value, ...) {
    counts(object, ...) <- value
    return(object)
})


#' @rdname getter_setter_functions
#' @export
setMethod("preprocessed", "Outrider2DataSet", 
        function(object, normalized=FALSE) {
    if(!("preprocessed" %in% assayNames(object)) || 
       modelParams(object, "preprocessing") == "none"){
        return(observed(object, normalized=normalized))
    }
            
    if(isTRUE(normalized)){
        object <- preprocess(object, normalized=TRUE)
    } 
    return(assay(object, "preprocessed"))
    
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("preprocessed", "Outrider2DataSet", function(object, value, ...) {
    assay(object, "preprocessed", ...) <- value
    validObject(object)
    return(object)
})


#' @rdname getter_setter_functions
#' @export
setMethod("expected", "Outrider2DataSet", function(object, ...) {
    if(!("normalizationFactors" %in% assayNames(object))){
        stop("Expected values have not yet been computed.")
    }
    return(normalizationFactors(object))
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("expected", "Outrider2DataSet", function(object, value, ...) {
    normalizationFactors(object) <- value
    validObject(object)
    return(object)
})

