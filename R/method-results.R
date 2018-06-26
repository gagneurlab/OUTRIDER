
compileResults <- function(object, padjCutoff, zScoreCutoff, round, all){
    features <- c("pValue", "padjust", "zScore", "l2fc")
    if(!is(object, 'OutriderDataSet')){
        stop('Please provide an OutriderDataSet')
    }
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
    if(isTRUE(round)){
        round <- 2
    }
    if(isFALSE(all)){
        object <- object[aberrant(object, padjCutoff=padjCutoff, 
                zScoreCutoff=zScoreCutoff, by="gene") > 0,]
        if(dim(object)[1]==0){
            warning('No significant events: use all=TRUE to print all counts.')
            return(data.table(geneID='a', sampleID='a', pValue=0.1, padjust=0.1,
                    zScore=1.1, l2fc=0.1, rawcounts=1, normcounts=1, mu=0.1,
                    disp=1.2, meanCorrected=23.3, AberrantBySample=0,
                    AberrantByGene=1, padj_rank=23.5)[0])
        }
    }
    
    DTfitparameters <- as.data.table(as.matrix(mcols(
            object, use.names=TRUE)[c('mu', 'disp')]), keep.rownames="geneID")
    
    countdataDF <- as.data.table(counts(object), keep.rownames="geneID")
    countNormDF <- as.data.table(counts(object, normalized=TRUE), 
            keep.rownames="geneID")
    aberrantDF <- as.data.table(keep.rownames="geneID", aberrant(
            object, padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff))
    
    featureDTs <- lapply(features, function(f){
        dt <- as.data.table(assays(object)[[f]], keep.rownames='geneID')
    })
    names(featureDTs) <- features
    
    featureDTs <- append(featureDTs, list(
            rawcounts=countdataDF, normcounts=countNormDF, aberrant=aberrantDF))
    
    featureMeltDTs <- lapply(names(featureDTs), function(n){
        featureDTs[[n]] %>% melt(id.vars="geneID", value.name=n,
                measure.vars=2:(dim(object)[2]+1), variable.name="sampleID")})
    
    tidyresults <- Reduce(x=featureMeltDTs, f=function(x, y){
        merge(x, y, by=c("geneID", "sampleID")) })
    
    tidyresults <- merge(tidyresults, DTfitparameters, by=c('geneID'))
    
    tidyresults <- merge(tidyresults, data.table(geneID=rownames(object), 
            meanCorrected=rowMeans(counts(object, normalized=TRUE))), 
            by=c('geneID'))
    
    tidyresults <- merge(tidyresults, data.table(sampleID=colnames(object), 
            AberrantBySample=aberrant(object, by='sample')), 
            by=c('sampleID'))
    
    tidyresults <- merge(tidyresults, data.table(geneID=rownames(object), 
            AberrantByGene=aberrant(object, by='gene')), 
            by=c('geneID'))
    
    if(isFALSE(all)){
        tidyresults <- tidyresults[aberrant == TRUE]
    }
    tidyresults[,padj_rank:=rank(padjust), by=sampleID]
    tidyresults <- tidyresults[order(padjust)]
    
    if(is.numeric(round)){
        col2round <- c("normcounts", "mu", "zScore", "l2fc", "disp",
                "meanCorrected")
        devNull <- sapply(col2round, function(x){
                tidyresults[,c(x):=round(get(x), as.integer(round))] })
    }
    return(tidyresults)
}

#' 
#' Accessor function for the 'results' object in an OutriderDataSet object. 
#' 
#' This function assembles a results table of significant outlier events based
#' on the given filter criteria. The table contains various information 
#' accumulated over the analysis pipeline. 
#' 
#' @param object An OutriderDataSet
#' @param padjCutoff The significant theshold to be applied
#' @param zScoreCutoff If provided additionally a z score threashold is applied
#' @param round Can be TRUE, defaults to 2, or an integer used for rounding
#'             with \code{\link[base]{round}} to make the output
#'             more user friendly
#' @param all By default FALSE, only significant read counts are listed in the 
#'             results. If TRUE all results are assembled resulting in a 
#'             data.table of length samples x genes
#' @param ... Additional arguments, currently not used
#' 
#' @return A data.table where each row is an outlier event and the columns
#'             contain additional information about this event. Eg padj, l2fc
#' 
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet()
#' ods <- OUTRIDER(ods)
#' res <- results(ods, all=TRUE)
#' 
#' res
#' 
#' @docType methods
#' @name results
#' @rdname results
#' @export results
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))

#' @rdname results
#' @export
setMethod("results", "OutriderDataSet", function(object, 
                    padjCutoff=0.05, zScoreCutoff=0, round=2, all=FALSE, 
                    ...){
    compileResults(object, padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff,
            round=round, all=all, ...)
})
