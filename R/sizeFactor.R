#'
#' SizeFactors accessor and estimation function
#' 
#' Accessor functions for the 'sizeFactors' information in a OutriderDataSet
#' object. 
#' 
#' The estimation of the size factors can make also use of existing 
#' log geometric means in the object. Those can be loaded from an 
#' existing model.
#'    
#' @docType methods
#' @name sizeFactors
#' @rdname sizeFactors
#' @aliases sizeFactors sizeFactors,OutriderDataSet-method 
#'         sizeFactors<- sizeFactors<-,OutriderDataSet,numeric-method
#'         estimateSizeFactors estimateSizeFactors,OutriderDataSet-method
#'         
#' @param object OutriderDataSet
#' @param value A numberic vector of sizeFactors
#' @return An OutriderDatasSet with the estimated sizeFactors or with the 
#'             getter function it returns a numeric vector containing the 
#'             sizeFactors.
#' 
#' @seealso \code{\link[DESeq2]{estimateSizeFactors}}
#' 
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' head(sizeFactors(ods))
#' 
#' sizeFactors(ods) <- runif(dim(ods)[2], 0.5, 1.5)
#' sizeFactors(ods)
#' counts(ods, normalized=TRUE)[1:10,1:10]
#'         
NULL

#' 
#' Getter function for the size factors
#' @noRd
sizeFactors.OutriderDataSet <- function(object){
    if (!"sizeFactor" %in% names(colData(object))){
        return(NULL)
    }
    sf <- colData(object)$sizeFactor
    names(sf) <- colnames(object)
    sf
}

#' @rdname sizeFactors
#' @export sizeFactors
setMethod("sizeFactors", signature(object="OutriderDataSet"), 
        sizeFactors.OutriderDataSet)


#'
#' Replace function for the size factors
#' @noRd
sizeFactors.replace.OutriderDataSet <- function(object, value){
    stopifnot(all(!is.na(value)))
    stopifnot(all(is.finite(value)))
    stopifnot(all(value > 0))
    
    colData(object)$sizeFactor <- value
    validObject( object )
    object
}

#' @rdname sizeFactors
#' @export "sizeFactors<-"
setReplaceMethod("sizeFactors", signature(object="OutriderDataSet", 
        value="numeric"), sizeFactors.replace.OutriderDataSet)


#'
#' OUTRIDER size Factor method which stores the log geom means.
#' @noRd
estimateSizeFactors.OUTRIDER <- function(object){
    if(!'loggeomeans' %in% names(mcols(object))){
        mcols(object)[['loggeomeans']] <- rowMeans(log(counts(object)))
    }
    loggeomeans <- mcols(object)[['loggeomeans']]
    
    if(all(is.infinite(loggeomeans))){
        stop(paste("Every gene contains at least one zero,",
                "cannot compute log geometric means"))
    }
    
    sf <- apply(counts(object), 2, function(cnts) {
        exp(median((log(cnts) - loggeomeans)[
            is.finite(loggeomeans) & cnts > 0]))
    })
    
    sizeFactors(object) <- sf
    validObject(object)
    return(object)
}

#' @rdname sizeFactors
#' @export estimateSizeFactors
setMethod("estimateSizeFactors", "OutriderDataSet", 
        estimateSizeFactors.OUTRIDER)

