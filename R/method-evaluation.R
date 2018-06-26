#' 
#' Number of aberrant events
#' 
#' Identifies the aberrant events and returns the number of aberrant counts per
#' gene or sample or returns a matrix indicating aberrant events.
#' 
#' @param ods An OutriderDataSet object
#' @param padjCutoff The padjust cutoff 
#' @param zScoreCutoff The absolute Z-score cutoff, 
#'             if NA or NULL no Z-score cutoff is used
#' @param by if the results should be summarized by 'sample', 
#'             'gene' or not at all (default).
#'
#' @return The number of aberrent events by gene or sample or a TRUE/FALSE 
#'       matrix of the size sample x gene of aberrent events.
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(n=50)
#' ods <- OUTRIDER(ods)
#' 
#' aberrant(ods)[1:10,1:10]
#' aberrant(ods, by="sample")
#' aberrant(ods, by="gene")
#' 
#' @rdname aberrant
#' 
#' @export
aberrant <- function(ods, padjCutoff=0.05, zScoreCutoff=0, 
                    by=c("none", "sample", "gene")){
    aberrantEvents <- assays(ods)[['padjust']] < padjCutoff
    if(isScalarNumeric(zScoreCutoff, na.ok=FALSE)){
        aberrantEvents <- aberrantEvents & 
                abs(assays(ods)[['zScore']]) > zScoreCutoff
    }
    
    return(switch(match.arg(by),
            none = aberrantEvents,
            sample = colSums(aberrantEvents),
            gene = rowSums(aberrantEvents)
    ))
}