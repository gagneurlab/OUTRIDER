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
setMethod("dispersions", "Outrider2DataSet",
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


observed.OUTRIDER2 <- function(object, ...){
    
    if(!("observed" %in% assayNames(object))){
        if(!("counts" %in% assayNames(object))){
            stop("Assay 'observed' does not exist.")
        } else{
            obs <- assay(object, "counts", ...)
        }
    } else{
        obs <- assay(object, "observed", ...)
    }
    
    return(obs)

}

#' @rdname getter_setter_functions
#' @export
setMethod("observed", "Outrider2DataSet", observed.OUTRIDER2)

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
        function(object, normalized=FALSE, minE=0.5) {
    if(!("preprocessed" %in% assayNames(object)) ){
        stop("Calculate preprocessed data first with ods <- preprocess(ods)")
    }
            
    prepro   <- assay(object, "preprocessed")
    if(!normalized){
        return(prepro)    
    }
            
    if(!is.null(normalizationFactors(object))){
        if(all(prepro >= 0, na.rm=TRUE) && 
                !("expectedMean" %in% colnames(mcols(object)))){
            if(!"expectedLogGeomMean" %in% colnames(mcols(object))){
                stop("The expectedLogGeomMean is missing in mcols(ods). ",
                    "Did you set the normaliationFactors by hand? ",
                    "Please use normalizationFactors(object) <- values.")
            }
            
            # use cached expected log geom mean values
            eMat <- pmax(normalizationFactors(object), minE)
            eLGM <- mcols(object)[["expectedLogGeomMean"]]
            return(prepro/eMat * eLGM)
        } else{
            # if negative values in matrix
            if(!"expectedMean" %in% colnames(mcols(object))){
                stop("The expectedMean is missing in mcols(ods). ",
                    "Did you set the normaliationFactors by hand? ",
                    "Please use normalizationFactors(object) <- values.")
            }
            eMat <- normalizationFactors(object)
            eLGM <- mcols(object)[["expectedMean"]] 
            return((prepro - eMat) + eLGM)
        }
    }
    
    # normalization by sizeFactors
    if(is.null(sizeFactors(object)) || any(is.na(sizeFactors(object)))) {
        stop(paste("first calculate size factors, add normalizationFactors,",
                    "or set normalized=FALSE"))
    }
    return(t(t(prepro) / sizeFactors(object)))
        
    
})

#' @rdname getter_setter_functions
#' @export
setReplaceMethod("preprocessed", "Outrider2DataSet", 
    function(object, value, ...) {
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

#' @rdname getter_setter_functions
#' @export
profile <- function(ods){
    return(slot(ods, "profile"))
}

#' @rdname getter_setter_functions
#' @export
`profile<-` <- function(ods, ..., value){
    stopifnot(value %in% c("outrider", "protrider", "other"))
    slot(ods, "profile") <- value
    return(ods)
}

#' Default paramenter getter functions
#' 
#' This is a collection of functions for easy access to the default parameter 
#' settings within the OUTRIDER2 model based on the profile specifed in the 
#' Outrider2DataSet. If needed, the user can adapt the parameters and pass 
#' them on to the OUTRIDER function.
#' 
#' @param ods An OutriderDataSet object.
#' @param ... Further arguments passed on to the underlying assay function.
#' @return A list of parameters. 
#' 
#' @name default_param_getters
#' @rdname default_param_getters
#' @aliases getDefaultPreproParams, getDefaultPvalueParams
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(10, 10)
#' getDefaultPreproParams()
#' getDefaultPvalueParams()
#' 
#' 
NULL

#' @rdname default_param_getters
#' @export
getDefaultPreproParams <- function(ods){
    
    profile <- profile(ods) # profile=c("outrider", "protrider", "other")
    params <- list()
    
    # preprocessing options
    params[["sf_norm"]] <- TRUE
    params[["noise_factor"]] <- 0.0 # only for python implementation
    
    if(profile == "outrider"){
        params[["data_trans"]] <- "log1p"
        params[["reverse_data_trans"]] <- "exp"
        params[["distribution"]] <- "nb"
    } else{
        params[["data_trans"]] <- NULL
        params[["reverse_data_trans"]] <- NULL
        params[["distribution"]] <- "gaussian"
    } 
    
    if(profile == "protrider"){
        params[["prepro_func"]] <- function(x) log2(x+1)
    } else{
        params[["prepro_func"]] <- NULL 
    }
    
    if(profile == "other"){
        params[["sf_norm"]] <- FALSE
        # for vst: prepro_func <- DESeq2::varianceStabilizingTransformation    
    }
    
    return(params)
}

#' @rdname default_param_getters
#' @export
getDefaultPvalueParams <- function(ods){
    
    profile <- profile(ods) # profile=c("outrider", "protrider", "other")
    
    params <- list()
    
    # options for pvalue calculation
    params[["fdr_method"]] <- "BY"
    params[["pvalue_alternative"]] <- "two.sided"
    
    if(profile == "outrider"){
        params[["distribution"]] <- "nb"
        params[["effect_types"]] <- c("fold_change", "zscores")
    } else if(profile == "protrider"){
        params[["distribution"]] <- "gaussian"
        params[["effect_types"]] <- c("fold_change", "zscores")
    } else{
        params[["distribution"]] <- "gaussian"
        params[["effect_types"]] <- c("zscores", "delta")
    }
    
    return(params)
}

