compileResults.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                    round=2, all=FALSE, ...){
    
    #
    # input check and parsing
    # 
    checkOutriderDataSet(object)
    checkFullAnalysis(object)
    
    if(is.null(rownames(object))){
        rownames(object) <- paste('feature', seq_len(nrow(object)), sep='_')
    }
    
    if(is.null(colnames(object))){
        colnames(object) <- paste('sample', seq_len(ncol(object)), sep='_')
    }
    
    if(isTRUE(round)){
        round <- 2
    }
    
    if(isFALSE(all)){
        abByGene <- aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="gene")
        abBySample <- aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="sample")
        object <- object[abByGene > 0, abBySample > 0]
    }
    
    if(nrow(object) == 0){
        if(isFALSE(all)){
            warning('No significant events: use all=TRUE to print all events.')
        } else {
            warning('Please provide an object with at least one feature.')
        }
        return(data.table(
                geneID=NA_character_,
                sampleID=NA_character_,
                pValue=NA_real_,
                padjust=NA_real_,
                zScore=NA_real_,
                l2fc=NA_real_,
                rawcounts=NA_integer_,
                normcounts=NA_real_,
                meanCorrected=NA_real_,
                theta=NA_real_,
                aberrant=NA,
                AberrantBySample=NA_integer_,
                AberrantByGene=NA_integer_,
                padj_rank=NA_real_)[0])
    }
    
    #
    # extract data
    # 
    ans <- data.table(
        geneID           = rownames(object), 
        sampleID         = rep(colnames(object), each=nrow(object)),
        pValue           = c(pValue(object)),
        padjust          = c(padj(object)),
        zScore           = c(zScore(object)),
        l2fc             = c(assay(object, "l2fc")),
        rawcounts        = c(counts(object)),
        normcounts       = c(counts(object, normalized=TRUE)),
        meanCorrected    = rowMeans(counts(object, normalized=TRUE)),
        theta            = theta(object),
        aberrant         = c(aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff)),
        AberrantBySample = rep(each=nrow(object), aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="sample")),
        AberrantByGene   = aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="gene"),
        padj_rank        = c(apply(padj(object), 2, rank)))
    
    # round columns if requested
    if(is.numeric(round)){
        devNull <- lapply(
                c("normcounts", "zScore", "l2fc", "theta", "meanCorrected"),
                function(x) ans[,c(x):=round(get(x), as.integer(round))] )
    }
    
    # 
    # keep only aberrent events and sort by padj value
    # 
    if(isFALSE(all)){
        ans <- ans[aberrant == TRUE]
    }
    ans <- ans[order(padjust)]
    
    return(ans)
}

#' 
#' Accessor function for the 'results' object in an OutriderDataSet object. 
#' 
#' This function assembles a results table of significant outlier events based
#' on the given filter criteria. The table contains various information 
#' accumulated over the analysis pipeline. 
#' 
#' @param object An OutriderDataSet object
#' @param padjCutoff The significant threshold to be applied
#' @param zScoreCutoff If provided additionally a z score threshold is applied
#' @param round Can be TRUE, defaults to 2, or an integer used for rounding
#'             with \code{\link[base]{round}} to make the output
#'             more user friendly
#' @param all By default FALSE, only significant read counts are listed in the 
#'             results. If TRUE all results are assembled resulting in a 
#'             data.table of length samples x genes
#' @param ... Additional arguments, currently not used
#' 
#' @return A data.table where each row is an outlier event and the columns
#'    contain additional information about this event. In details the 
#'    table contains: 
#'    \itemize{
#'      \item{geneID: }{The gene ID as provided by the user.}
#'      \item{sampleID: }{The sample ID as provided by the user.}
#'      \item{pValue: }{The nominal P-value as computed by 
#'            \code{\link[OUTRIDER]{computePvalues}}.}
#'      \item{padjust: }{The FDR corrected P-value as computed by 
#'            \code{\link[OUTRIDER]{computePvalues}}.}
#'      \item{zScore: }{The z score as computed by 
#'            \code{\link[OUTRIDER]{computeZscores}}.}
#'      \item{l2fc: }{The log\eqn{_2}{[2]} fold change as computed by 
#'            \code{\link[OUTRIDER]{computeZscores}}.}
#'      \item{rawcounts: }{The observed read counts.}
#'      \item{normcounts: }{The expected count given the fitted 
#'            autoencoder model.}
#'      \item{meanCorrected: }{For this gene, the mean of the expected counts
#'           given the fitted autoencoder model.}
#'      \item{theta: }{The dispersion parameter of the NB distribution 
#'            for the given gene.}
#'      \item{aberrant: }{The outlier status of this event: 
#'            \code{TRUE} or \code{FALSE}.}
#'      \item{AberrantBySample: }{Number of outliers for the given sample.}
#'      \item{AberrantByGene: }{Number of outliers for the given gene.}
#'      \item{padj_rank: }{Rank of this outlier within the given sample.}
#'    }
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
#' @name results
#' @rdname results
#' @aliases results results,OutriderDataSet-method
#' 
#' @export
setMethod("results", "OutriderDataSet", compileResults.OUTRIDER)

