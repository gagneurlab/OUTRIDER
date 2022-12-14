#' Getter/Setter functions
#' 
#' This is a collection of small accessor/setter functions for easy access to
#' the values within the OUTRIDER model.
#' 
#' @param ods,object An OutriderDataSet object.
#' @param subsetName Name of a gene subset for which to store ore retrieve 
#'              FDR corrected p values 
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
padj <- function(ods, subsetName=NULL){
    aname <- 'padjust'
    if(!is.null(subsetName)){
        aname <- paste0(aname, "_", subsetName)
    }
    if(!aname %in% assayNames(ods)){
        stop('Please compute first the P-values before retrieving', '
                the adjusted ones.')
    }
    assay(ods, aname)
}

`padj<-` <- function(ods, subsetName=NULL, ..., value){
    stopifnot(is.matrix(value))
    stopifnot(dim(ods) == dim(value))
    aname <- 'padjust'
    if(!is.null(subsetName)){
        aname <- paste0(aname, "_", subsetName)
    }
    assay(ods, aname, ...) <- value
    return(ods)
}


#' @rdname getter_setter_functions
#' @aliases dispersions dispersions,OutriderDataSet-method
#' @seealso \code{\link[DESeq2]{estimateDispersions}}
#' @export
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
