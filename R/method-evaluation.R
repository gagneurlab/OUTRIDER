aberrant.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                    by=c("none", "sample", "gene"), genesToTest=NULL, ...){
    checkFullAnalysis(object)
    
    dots <- list(...)
    if(is.null(genesToTest) && is.null(dots[["FDR_subset"]])){
        aberrantEvents <- padj(object) <= padjCutoff
    } else{
        if("FDR_subset" %in% names(dots)){
            fdr_subset <- dots[["FDR_subset"]]
        } else{
            subsetName <- "subset"
            object <- padjOnSubset(ods=object, genesToTest=genesToTest, 
                                    subsetName=subsetName)
            fdr_subset <- metadata(object)[[paste("FDR", subsetName, sep="_")]]
        }
        
        # define aberrant status based on whether genes/sample tuples are 
        # part of the given subset
        aberrantEvents <- matrix(FALSE, nrow=nrow(object), ncol=ncol(object))
        rownames(aberrantEvents) <- rownames(object)
        colnames(aberrantEvents) <- colnames(object)
        subset_idx <- as.matrix(
            fdr_subset[sampleID %in% colnames(object) & 
                        FDR_subset <= padjCutoff, 
                            .(idx, sapply(sampleID, 
                                function(s) which(colnames(object) == s)))]
        )
        if(nrow(subset_idx) > 0){
            aberrantEvents[subset_idx] <- TRUE
        }
    }
    
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
#' @param genesToTest A named list giving a subset of genes per sample to which  
#'              FDR correction should be restricted before determining aberrant 
#'              status. The names of the list must correspond to the sampleIDs 
#'              in the ods object. See examples for how to use this feature.
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
#' # example of retrieving aberrant status with FDR correction limited to a 
#' # set of genes of interest
#' genesToTest <- list("sample_1"=sample(rownames(ods)[1:10], 3), 
#'                     "sample_2"=sample(rownames(ods)[1:10], 8), 
#'                     "sample_6"=sample(rownames(ods)[1:10], 5))
#' aberrant(ods, genesToTest=genesToTest)[1:10,1:10]
#' tail(sort(aberrant(ods, genesToTest=genesToTest, by="sample")))
#' tail(sort(aberrant(ods, genesToTest=genesToTest, by="gene")))
#' 
#' @rdname aberrant
#' @export
setMethod("aberrant", signature="OutriderDataSet", aberrant.OUTRIDER)
