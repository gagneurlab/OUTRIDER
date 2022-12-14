aberrant.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                    by=c("none", "sample", "gene"), subsetName=NULL, ...){
    checkFullAnalysis(object)
    
    aberrantEvents <- padj(object, subsetName=subsetName) <= padjCutoff
    
    if(isScalarNumeric(zScoreCutoff, na.ok=FALSE)){
        aberrantEvents <- aberrantEvents & abs(zScore(object)) >= zScoreCutoff
    }
    
    return(switch(match.arg(by),
            none = aberrantEvents,
            sample = colSums(aberrantEvents, na.rm=TRUE),
            gene = rowSums(aberrantEvents, na.rm=TRUE)
    ))
}

#' 
#' Number of aberrant events
#' 
#' Identifies the aberrant events and returns the number of aberrant counts per
#' gene or sample or returns a matrix indicating aberrant events.
#' 
#' @param object An OutriderDataSet object
#' @param padjCutoff The padjust cutoff 
#' @param zScoreCutoff The absolute Z-score cutoff, 
#'             if NA or NULL no Z-score cutoff is used
#' @param by if the results should be summarized by 'sample', 
#'             'gene' or not at all (default).
#' @param subsetName The name of a subset of genes for which FDR values have 
#'             been previously computed. Those FDR values on the subset will 
#'             then be used to determine aberrant status.
#' @param ... Currently not in use.
#'
#' @return The number of aberrent events by gene or sample or a TRUE/FALSE 
#'             matrix of the size sample x gene of aberrent events.
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- OUTRIDER(ods, implementation='pca')
#' 
#' aberrant(ods)[1:10,1:10]
#' tail(sort(aberrant(ods, by="sample")))
#' tail(sort(aberrant(ods, by="gene")))
#' 
#' @rdname aberrant
#' @export
setMethod("aberrant", signature="OutriderDataSet", aberrant.OUTRIDER)
