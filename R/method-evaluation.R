## Evaluation functions


#' aberrant
#' 
#' identifies the aberrant counts or returns the number of aberrant counts per
#' gene or sample.
#' 
#' @param ods a OutriderDataSet object
#' @param padj the padjusted cutoff 
#' @param zScore the absolute Z-score cutoff
#' @param by if the results should be summarized by 'sample', 'gene' or not 
#' at all.
#'
#' @return the number of aberrent events by gene or sample or a TRUE/FALSE 
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
aberrant <- function(ods, padj=0.05, zScore=3, 
                     by=c("none", "sample", "gene")){
    
    aberrantEvents <- assays(ods)[['padjust']] < padj &
        abs(assays(ods)[['zScore']]) > zScore
    
    return(switch(match.arg(by),
                  none = aberrantEvents,
                  sample = colSums(aberrantEvents),
                  gene = rowSums(aberrantEvents)
    ))
}


significant <- function(ods, padj=0.05, 
                        by=c("none", "sample", "gene")){
    
    signifEvents <- assays(ods)[['padjust']] < padj
    
    return(switch(match.arg(by),
                  none = signifEvents,
                  sample = colSums(signifEvents),
                  gene = rowSums(signifEvents)
    ))
}
