# compileResults.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
#                     round=2, all=FALSE, ...){
#     
#     #
#     # input check and parsing
#     # 
#     checkOutriderDataSet(object)
#     checkFullAnalysis(object)
#     
#     if(is.null(rownames(object))){
#         rownames(object) <- paste('feature', seq_len(nrow(object)), sep='_')
#     }
#     
#     if(is.null(colnames(object))){
#         colnames(object) <- paste('sample', seq_len(ncol(object)), sep='_')
#     }
#     
#     if(isTRUE(round)){
#         round <- 2
#     }
#     
#    if(isFALSE(all)){
#        abByGene <- aberrant(object,
#                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff,
#                by="feature")
#        abBySample <- aberrant(object,
#                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="sample")
#        object <- object[abByGene > 0, abBySample > 0]
#    }
#     
#     if(nrow(object) == 0){
#         if(isFALSE(all)){
#             warning('No significant events: use all=TRUE to print all events.')
#         } else {
#             warning('Please provide an object with at least one feature.')
#         }
#         return(data.table(
#                 geneID=NA_character_,
#                 sampleID=NA_character_,
#                 pValue=NA_real_,
#                 padjust=NA_real_,
#                 zScore=NA_real_,
#                 l2fc=NA_real_,
#                 rawcounts=NA_integer_,
#                 normcounts=NA_real_,
#                 meanCorrected=NA_real_,
#                 theta=NA_real_,
#                 aberrant=NA,
#                 AberrantBySample=NA_integer_,
#                 AberrantByGene=NA_integer_,
#                 padj_rank=NA_real_)[0])
#     }
#     
#     #
#     # extract data
#     # 
#     ans <- data.table(
#         geneID           = rownames(object), 
#         sampleID         = rep(colnames(object), each=nrow(object)),
#         pValue           = c(pValue(object)),
#         padjust          = c(padj(object)),
#         zScore           = c(zScore(object)),
#         l2fc             = c(assay(object, "l2fc")),
#         rawcounts        = c(counts(object)),
#         normcounts       = c(counts(object, normalized=TRUE)),
#         meanCorrected    = rowMeans(counts(object, normalized=TRUE)),
#         theta            = variability(object), # return theta for NB
#         aberrant         = c(aberrant(object, 
#                 padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff)),
#         AberrantBySample = rep(each=nrow(object), aberrant(object, 
#                   padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, 
#                   by="sample")),
#         AberrantByGene   = aberrant(object, 
#                   padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, 
#                   by="feature"),
#         padj_rank        = c(apply(padj(object), 2, rank)))
#     
#     # round columns if requested
#     if(is.numeric(round)){
#         devNull <- lapply(
#                 c("normcounts", "zScore", "l2fc", "theta", "meanCorrected"),
#                 function(x) ans[,c(x):=round(get(x), as.integer(round))] )
#     }
#     
#     # drop theta column if not NB used
#     if(modelParams(object, "distribution") != "negative binomial"){
#         ans[,theta:=NULL]
#     }
#     
#     # 
#     # keep only aberrent events and sort by padj value
#     # 
#     if(isFALSE(all)){
#         ans <- ans[aberrant == TRUE]
#     }
#     ans <- ans[order(padjust)]
#     
#     return(ans)
# }

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
#' @name results
#' @rdname results
#' @aliases results results,OutriderDataSet-method
#' 
#' @export
# setMethod("results", "OutriderDataSet", compileResults.OUTRIDER)


#' internal result function
#' @noRd
compileResults.OUTRIDER2 <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                                    round=2, all=FALSE, ...){
    
    #
    # input check and parsing
    # 
    checkOutrider2DataSet(object)
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
        abByFeature <- aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="feature")
        abBySample <- aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="sample")
        object <- object[abByFeature > 0, abBySample > 0]
    }
    
    if(nrow(object) == 0){
        if(isFALSE(all)){
            warning('No significant events: use all=TRUE to print all events.')
        } else {
            warning('Please provide an object with at least one feature.')
        }
        return(data.table(
            featureID=NA_character_,
            sampleID=NA_character_,
            pValue=NA_real_,
            padjust=NA_real_,
            zScore=NA_real_,
            delta=NA_real_,
            observed=NA_real_,
            preprocessed=NA_real_,
            # expected=NA_real_,
            normalized=NA_real_,
            meanCorrected=NA_real_,
            sd=NA_real_,
            pvalDistribution=NA_real_,
            aberrant=NA,
            AberrantBySample=NA_integer_,
            AberrantByFeature=NA_integer_,
            padj_rank=NA_real_)[0])
    }
    
    #
    # extract data
    # 
    ans <- data.table(
        featureID        = rownames(object), 
        sampleID         = rep(colnames(object), each=nrow(object)),
        pValue           = c(pValue(object)),
        padjust          = c(padj(object)),
        zScore           = c(zScore(object)),
        effect           = c(effect(object)),
        observed         = c(observed(object)),
        preprocessed     = c(preprocessed(object)),
        # expected         = c(expected(object)),
        normalized       = c(observed(object, normalized=TRUE)),
        meanCorrected    = rowMeans(observed(object, normalized=TRUE)),
        sd               = variability(object),
        pvalDistribution = modelParams(object, "distribution"),
        aberrant         = c(aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff)),
        AberrantBySample = rep(each=nrow(object), aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="sample")),
        AberrantByFeature   = aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="feature"),
        padj_rank        = c(apply(padj(object), 2, rank)))
    
    # round columns if requested
    if(is.numeric(round)){
        devNull <- lapply(
            c("normalized", "zScore", "effect", "sd", "meanCorrected", 
              "observed", "preprocessed"),
            function(x) ans[,c(x):=round(get(x), as.integer(round))] )
    }
    
    # rename columns for neg bin distribution
    if(modelParams(object, "distribution") == "negative binomial"){
        setnames(ans, "sd", "theta")
        setnames(ans, "observed", "rawcounts")
        setnames(ans, "normalized", "normcounts")
        setnames(ans, "effect", "l2fc")
    }
    # rename effect column to delta for gaussian
    if(modelParams(object, "distribution") == "gaussian"){
        setnames(ans, "effect", "delta")
    }
    
    # remove preprocess column if no preprocesssing was done
    if(modelParams(object, "preprocessing") == "none"){
        ans[,preprocessed:=NULL]
    }
    
    # set featureID -> geneID for OutriderDataSet
    if(is(object, "OutriderDataSet")){
        setnames(ans, "featureID", "geneID")
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

#' @name results
#' @rdname results
#' @aliases results results,Outrider2DataSet-method
#' 
#' @export
setMethod("results", "Outrider2DataSet", compileResults.OUTRIDER2)
