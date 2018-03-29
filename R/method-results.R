
compileResults <- function(object, padj=0.05, zScore=3, 
                    round=TRUE, all=FALSE){
    padjCut <- padj
    zScoreCut <- zScore
    features <- c("pValue", "padjust", "zScore", "l2fc")
    if(!"pValue" %in% assayNames(object)){
        stop(paste0("The P-values are not computed yet. Please run the ",
                "following command:\n\tods <- computePvalues(ods)"))
    }
    if(!"zScore" %in% assayNames(object)){
        stop(paste0("The Z-scores are not computed yet. Please run the ",
                "following command:\n\tods <- computeZscores(ods)"))
    }
    if(!all(features %in% assayNames(object))){
        stop(paste0("Some of the assays are missing. Please run the full ", 
                "analysis before extracting the results.",
                "\n\tods <- OUTRIDER(ods)"))
    }
    
    if(all==FALSE){
        object <- object[rowMins(assays(object)[['padjust']]) <= padjCut]
        if(dim(object)[1]==0){
            warning('No significant events: use all=TRUE to print all counts.')
        }
    }
    
    DTfitparameters <- as.data.table(as.matrix(mcols(
            object, use.names=TRUE)[c('mu', 'disp')]), keep.rownames="geneID")
    
    countdataDF <- as.data.table(counts(object), keep.rownames="geneID")
    countNormDF <- as.data.table(counts(object, normalized=TRUE), 
            keep.rownames="geneID")
    
    featureDTs <- lapply(features, function(f){
        dt <- as.data.table(assays(object)[[f]], keep.rownames='geneID')
    })
    names(featureDTs) <- features
    
    featureDTs <- append(featureDTs, list(
            rawcounts=countdataDF, normcounts=countNormDF))
    
    featureMeltDTs <- lapply(names(featureDTs), function(n){
        featureDTs[[n]] %>% melt(id.vars="geneID", value.name=n,
                measure.vars=2:(dim(object)[2]+1), variable.name="sampleID")})
    
    tidyresults <- featureMeltDTs[[1]]
    for(i in 2:length(featureMeltDTs)){
        tidyresults <- merge(tidyresults, featureMeltDTs[[i]], 
                by=c("geneID", "sampleID"))
    }
    
    tidyresults <- merge(tidyresults, DTfitparameters, by=c('geneID'))
    
    tidyresults <- merge(tidyresults, data.table(geneID = rownames(object), 
            meanCorrected = rowMeans(counts(object, normalized=TRUE))), 
            by=c('geneID'))
    
    tidyresults <- merge(tidyresults, data.table(sampleID=colnames(object), 
            AberrantBySample = aberrant(object, by='sample')), 
            by=c('sampleID'))
    
    tidyresults <- merge(tidyresults, data.table(geneID = rownames(object), 
            AberrantByGene = aberrant(object, by='gene')), 
            by=c('geneID'))
    
    if(all==FALSE){
        tidyresults <- tidyresults[
                padjust <= padjCut & abs(zScore) >= zScoreCut][order(padjust)]
    }
    tidyresults[,padj_rank := rank(padjust), by = sampleID]
    
    if(round == TRUE){
        tidyresults[,c("normcounts", "mu", "zScore", "l2fc", "disp"):=list(
                round(normcounts, 2), round(mu, 2), round(zScore,2), 
                round(l2fc, 2), round(disp, 2))]
    }
    return(tidyresults)
}

#' 
#' Accessor functions for the 'result' object in a OutriderDataSet
#' object. 
#' 
#' This function extracts all results based on the given filter criteria.
#' 
#' @param object OutriderDataSet
#' @param padj padj cutoff
#' @param zScore absolute Z-score cutoff 
#' @param round rounding the results.
#' @param all By default FALSE, only significant read counts are listed in the 
#' results. 
#' @param ... additional arguments, currently not used
#' 
#' @return A data.table where each row is an outlier event and the columns
#'             contain additional information about this event. Eg padj, l2fc
#' 
#' @examples
#' set.seed(42)
#' ods <- makeExampleOutriderDataSet()
#' ods <- OUTRIDER(ods)
#' res <- results(ods)
#' 
#' res
#' 
#' @docType methods
#' @name results
#' @rdname results
#' @aliases results results,OutriderDataSet-method
#' @export results
setMethod("results", "OutriderDataSet", compileResults)
