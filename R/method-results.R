
compileResults <- function(object, padjCutoff, zScoreCutoff, round, all, 
                           BPPARAM=bpparam()){
    checkOutriderDataSet(object)
    checkFullAnalysis(object)
    
    features <- c("pValue", "padjust", "zScore", "l2fc")
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
            return(data.table(geneID=NA_character_, sampleID=NA_character_, 
                    pValue=NA_real_, padjust=NA_real_, zScore=NA_real_,
                    l2fc=NA_real_, rawcounts=NA_integer_, normcounts=NA_real_,
                    theta=NA_real_, meanCorrected=NA_real_, 
                    AberrantBySample=NA_integer_, AberrantByGene=NA_integer_,
                    padj_rank=NA_real_)[0])
        }
    }
    
    if(is.null(rownames(object))){
        rownames(object) <- paste('feature', seq_len(nrow(object)), sep='_')
    }
    DTfitparameters <- data.table(geneID=rownames(object), theta=theta(object))
    
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
    
    featureMeltDTs <- bplapply(names(featureDTs), dt=featureDTs, 
            BPPARAM=BPPARAM, mv=seq_len(ncol(object)) + 1,
            FUN=function(n, dt, mv){
                dt[[n]] %>% melt(id.vars="geneID", value.name=n,
                        measure.vars=mv, variable.name="sampleID")})
    
    tidyresults <- Reduce(x=featureMeltDTs, f=function(x, y){
        merge(x, y, by=c("geneID", "sampleID")) })
    
    # merge by geneID
    tidyresults <- merge(tidyresults, DTfitparameters, by=c('geneID'))
    tidyresults <- merge(tidyresults, data.table(geneID=rownames(object), 
            meanCorrected=rowMeans(counts(object, normalized=TRUE))), 
            AberrantByGene=aberrant(object, by='gene'), 
            by=c('geneID'))
    
    # merge by sampleID
    tidyresults <- merge(tidyresults, data.table(sampleID=colnames(object), 
            AberrantBySample=aberrant(object, by='sample')), 
            by=c('sampleID'))
    
    if(isFALSE(all)){
        tidyresults <- tidyresults[aberrant == TRUE]
    }
    tidyresults[,padj_rank:=rank(padjust), by=sampleID]
    tidyresults <- tidyresults[order(padjust)]
    
    if(is.numeric(round)){
        col2round <- c("normcounts", "zScore", "l2fc", "theta",
                "meanCorrected")
        devNull <- lapply(col2round, function(x){
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
#' \dontshow{
#'     ods <- ods[1:10,1:10]
#' }
#' ods <- OUTRIDER(ods)
#' 
#' res <- results(ods, all=TRUE)
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
