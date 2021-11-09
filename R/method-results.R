#' internal result function
#' @noRd
compileResults.OUTRIDER2 <- function(object, padjCutoff=0.05, 
                                    zScoreCutoff=NULL, l2fcCutoff=NULL, 
                                    deltaCutoff=NULL, round=2, all=FALSE, ...){
    
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
        abByFeature <- aberrant(object, padjCutoff=padjCutoff, 
                zScoreCutoff=zScoreCutoff, l2fcCutoff=l2fcCutoff, 
                deltaCutoff=deltaCutoff, by="feature")
        abBySample <- aberrant(object, padjCutoff=padjCutoff, 
                zScoreCutoff=zScoreCutoff, l2fcCutoff=l2fcCutoff, 
                deltaCutoff=deltaCutoff, by="sample")
        object <- object[abByFeature > 0, abBySample > 0]
    }
    
    if(nrow(object) == 0){
        if(isFALSE(all)){
            warning('No significant events: use all=TRUE to print all events.')
        } else {
            warning('Please provide an object with at least one feature.')
        }
        if(profile(ods) == "outrider"){
            return(data.table(
                geneID=NA_character_,
                sampleID=NA_character_,
                pValue=NA_real_,
                padjust=NA_real_,
                zScore=NA_real_,
                fc=NA_real_,
                log2fc=NA_real_,
                rawcounts=NA_real_,
                expected_counts=NA_real_,
                normcounts=NA_real_,
                meanCorrected=NA_real_,
                theta=NA_real_,
                sizefactor=NA_real_,
                pvalDistribution=NA_character_,
                aberrant=NA,
                AberrantBySample=NA_integer_,
                AberrantByFeature=NA_integer_,
                padj_rank=NA_real_)[0])
        } else{
            return(data.table(
                featureID=NA_character_,
                sampleID=NA_character_,
                pValue=NA_real_,
                padjust=NA_real_,
                zScore=NA_real_,
                fc=NA_real_,
                log2fc=NA_real_,
                delta=NA_real_,
                input_value=NA_real_,
                preprocessed_raw=NA_real_,
                preprocessed_expected=NA_real_,
                normalized=NA_real_,
                meanCorrected=NA_real_,
                sd=NA_real_, # theta
                sizefactor=NA_real_,
                pvalDistribution=NA_character_,
                aberrant=NA,
                AberrantBySample=NA_integer_,
                AberrantByFeature=NA_integer_,
                padj_rank=NA_real_)[0])
        }
        
    }
    
    #
    # extract data
    # 
    if("zScore" %in% assayNames(object)){
        zscores <- zScore(object)
    } else{
        zscores <- NA
    }
    if("l2fc" %in% assayNames(object)){
        l2fc <- assay(object, "l2fc")
    } else{
        l2fc <- NA
    }
    if("delta" %in% assayNames(object)){
        delta <- assay(object, "delta")
    } else{
        delta <- NA
    }
    
    ans <- data.table(
        featureID        = rownames(object), 
        sampleID         = rep(colnames(object), each=nrow(object)),
        pValue           = c(pValue(object)),
        padjust          = c(padj(object)),
        zScore           = c(zscores),
        fc               = c(2^l2fc),
        log2fc           = c(l2fc),
        delta            = c(delta), 
        input_value      = c(observed(object)),
        preprocessed_raw = c(preprocessed(object)),
        preprocessed_expected = c(expected(object)),
        normalized       = c(preprocessed(object, normalized=TRUE)),
        meanCorrected    = rowMeans(preprocessed(object, normalized=TRUE), 
                                    na.rm=TRUE),
        sd               = variability(object),
        sizefactor       = rep(each=nrow(object), sizeFactors(object)),
        pvalDistribution = metadata(object)$pvalDistribution,
        aberrant         = c(aberrant(object, padjCutoff=padjCutoff, 
                            zScoreCutoff=zScoreCutoff, l2fcCutoff=l2fcCutoff,
                            deltaCutoff=deltaCutoff)),
        AberrantBySample = rep(each=nrow(object), aberrant(object, 
                            padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, 
                            l2fcCutoff=l2fcCutoff, deltaCutoff=deltaCutoff, 
                            by="sample")),
        AberrantByFeature = aberrant(object, padjCutoff=padjCutoff, 
                            zScoreCutoff=zScoreCutoff, l2fcCutoff=l2fcCutoff,
                            deltaCutoff=deltaCutoff, by="feature"),
        padj_rank        = c(apply(padj(object), 2, rank)))
    
    # round columns if requested
    if(is.numeric(round)){
        devNull <- lapply(
            c("normalized", "zScore", "fc", "log2fc", "delta", "sd", 
                "meanCorrected", "input_value", "preprocessed_raw", 
                "preprocessed_expected", "sizefactor"),
            function(x) ans[,c(x):=round(get(x), as.integer(round))] )
    }
    
    # drop effect columns containing only NAs
    devNull <- lapply(c("zScore", "fc", "log2fc", "delta"), 
                    function(x){
                        vals <- ans[,unique(get(x))]
                        if(all(is.na(vals))){
                            ans[,(x):=NULL]
                        }
                    })
    
    # remove preprocess column if no preprocesssing was done
    if(is.null(metadata(object)$prepro_options$prepro_func)){
        ans[,preprocessed_raw:=NULL]
        setnames(ans, "preprocessed_expected", "expected_value")
    }
    
    # drop sizefactor column if sf_norm=FALSE (all equal to 1)
    if(all(ans[,sizefactor] == 1)){
        ans[,sizefactor:=NULL]
    }
    
    # rename columns for neg bin distribution
    if(metadata(object)$pvalDistribution == "nb"){
        if("preprocessed_raw" %in% colnames(ans)){
            setnames(ans, "preprocessed_raw", "rawcounts")
            setnames(ans, "preprocessed_expected", "expected_counts")
        } else{
            setnames(ans, "input_value", "rawcounts")
            setnames(ans, "expected_value", "expected_counts")
        }
        setnames(ans, "normalized", "normcounts")
        setnames(ans, "sd", "theta")
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
#' @param l2fcCutoff If provided additionally a threshold on the log2 fold 
#'             change is applied
#' @param deltaCutoff If provided additionally a threshold on delta is applied
#' @param round Either FALSE or an integer used for rounding 
#'             with \code{\link[base]{round}} to make the output more user 
#'             friendly. Defaults to 2.
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
#' @aliases results results,Outrider2DataSet-method
#' 
#' @export
setMethod("results", "Outrider2DataSet", compileResults.OUTRIDER2)
# #' @name results
# #' @rdname results
# #' @aliases results results,OutriderDataSet-method
# #' 
# #' @export
# setMethod("results", "OutriderDataSet", compileResults.OUTRIDER)

