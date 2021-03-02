# aberrant.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
#                     by=c("none", "sample", "gene")){
#     checkFullAnalysis(object)
#     
#     aberrantEvents <- padj(object) <= padjCutoff
#     if(isScalarNumeric(zScoreCutoff, na.ok=FALSE)){
#         aberrantEvents <- aberrantEvents & abs(zScore(object)) >= zScoreCutoff
#     }
#     
#     return(switch(match.arg(by),
#             none = aberrantEvents,
#             sample = colSums(aberrantEvents, na.rm=TRUE),
#             gene = rowSums(aberrantEvents, na.rm=TRUE)
#     ))
# }

aberrant.OUTRIDER2 <- function(object, padjCutoff=0.05, zScoreCutoff=NULL, 
                                l2fcCutoff=NULL, deltaCutoff=NULL,
                                by=c("none", "sample", "feature", "gene")){
    checkFullAnalysis(object)
    
    aberrantEvents <- padj(object) <= padjCutoff
    if(!is.null(zScoreCutoff) && isScalarNumeric(zScoreCutoff, na.ok=FALSE)){
        aberrantEvents <- aberrantEvents & abs(zScore(object)) >= zScoreCutoff
    }
    if(!is.null(l2fcCutoff) && isScalarNumeric(l2fcCutoff, na.ok=FALSE)){
        if(!("l2fc" %in% assayNames(object))){
            stop("Cannot apply l2fcCutoff because no l2fc values have been ", 
                "computed. Please add them first with ",
                "ods <- computeEffectSizes(ods, effect_types='fold_change')")
        }
        l2fc <- assay(object, "l2fc")
        aberrantEvents <- aberrantEvents & abs(l2fc) >= l2fcCutoff
    }
    if(!is.null(deltaCutoff) && isScalarNumeric(deltaCutoff, na.ok=FALSE)){
        if(!("delta" %in% assayNames(object))){
            stop("Cannot apply deltaCutoff because no delta values have been ", 
                "computed. Please add them first with ",
                "ods <- computeEffectSizes(ods, effect_types='delta')")
        }
        delta <- assay(object, "delta")
        aberrantEvents <- aberrantEvents & abs(delta) >= deltaCutoff
    }
    
    return(switch(match.arg(by),
                none = aberrantEvents,
                sample = colSums(aberrantEvents, na.rm=TRUE),
                feature = rowSums(aberrantEvents, na.rm=TRUE),
                gene = rowSums(aberrantEvents, na.rm=TRUE) # same as feature
    ))
}

#' 
#' Number of aberrant events
#' 
#' Identifies the aberrant events and returns the number of aberrant counts per
#' feature (e.g. gene) or sample or returns a matrix indicating aberrant events.
#' 
#' @param object An OutriderDataSet object
#' @param padjCutoff The padjust cutoff 
#' @param zScoreCutoff The absolute Z-score cutoff, 
#'             if NA or NULL no Z-score cutoff is used
#' @param by if the results should be summarized by 'sample', 
#'             'feature' (can also be called 'gene') or not at all (default).
#' @param ... Currently not in use.
#'
#' @return The number of aberrant events by feature (e.g. gene) or sample or 
#'             a TRUE/FALSE matrix of the size sample x gene of aberrant events.
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- OUTRIDER(ods, implementation='pca')
#' 
#' aberrant(ods)[1:10,1:10]
#' tail(sort(aberrant(ods, by="sample")))
#' tail(sort(aberrant(ods, by="feature")))
#' 
#' @rdname aberrant
#' @export
setMethod("aberrant", signature="Outrider2DataSet", aberrant.OUTRIDER2)
# setMethod("aberrant", signature="OutriderDataSet", aberrant.OUTRIDER)

