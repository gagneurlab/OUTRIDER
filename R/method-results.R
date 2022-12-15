compileResults.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                    round=2, all=FALSE, genesToTest=NULL, subsetName=NULL, 
                    ...){
    
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
    
    if(!is.null(genesToTest)){
        subsetName <- ifelse(is.null(subsetName), "subset", subsetName)
        FDR_set <- subsetName
        object <- padjOnSubset(ods=object, genesToTest=genesToTest, 
                               subsetName=subsetName)
    } else{
        FDR_set <- "transcriptome-wide"
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
                padj_rank=NA_real_,
                FDR_set=NA_character_)[0])
    }
    
    #
    # extract data
    #
    ans <- data.table(
        geneID           = rownames(object), 
        sampleID         = rep(colnames(object), each=nrow(object)),
        pValue           = c(pValue(object)),
        padjust          = c(padj(object, subsetName=subsetName)),
        zScore           = c(zScore(object)),
        l2fc             = c(assay(object, "l2fc")),
        rawcounts        = c(counts(object)),
        normcounts       = c(counts(object, normalized=TRUE)),
        meanCorrected    = rowMeans(counts(object, normalized=TRUE)),
        theta            = theta(object),
        aberrant         = c(aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff,
                subsetName=subsetName)),
        AberrantBySample = rep(each=nrow(object), aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="sample",
                subsetName=subsetName)),
        AberrantByGene   = aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="gene",
                subsetName=subsetName),
        padj_rank        = c(apply(padj(object), 2, rank)),
        FDR_set          = FDR_set)
    
    # round columns if requested
    if(is.numeric(round)){
        devNull <- lapply(
                c("normcounts", "zScore", "l2fc", "theta", "meanCorrected"),
                function(x) ans[,c(x):=round(get(x), as.integer(round))] )
    }
    
    # 
    # keep only aberrant events and sort by padj value
    # 
    if(isFALSE(all)){
        ans <- ans[aberrant == TRUE]
        
    } else if(!is.null(genesToTest)){
        # if return full subset is requested, retrieve those as all non-NA padj
        ans <- ans[!is.na(padjust),]
    }
    ans <- ans[order(padjust)]
    
    return(ans)
}

compileResultsAll.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                                    round=2, all=FALSE, subsets=NULL, 
                                    returnTranscriptomewideResults=TRUE, ...){
    if(is.null(subsets)){
        # always return transcriptome-wide results if no subsets are provided
        returnTranscriptomewideResults <- TRUE
    }
    
    # first retrieve transcriptome-wide results
    if(isTRUE(returnTranscriptomewideResults)){
        res <- compileResults.OUTRIDER(object=object, padjCutoff=padjCutoff, 
                            zScoreCutoff=zScoreCutoff, all=all, round=round,
                            genesToTest=NULL, subsetName=NULL, ...)
    }
    
    # add results for FDR_subsets if requested
    if(!is.null(subsets)){
        stopifnot(is.list(subsets))
        stopifnot(!is.null(names(subsets)))
        resls_subsets <- lapply(names(subsets), function(setName){
            geneList_sub <- subsets[[setName]]
            res_sub <- compileResults.OUTRIDER(object=object, 
                            padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, 
                            all=all, round=round, genesToTest=geneList_sub, 
                            subsetName=setName, ...)
            return(res_sub)
        })
        if(isTRUE(returnTranscriptomewideResults)){
            res <- rbindlist(append(list(res), resls_subsets))
        } else{
            res <- rbindlist(resls_subsets)
        }
        
        # sort it if existing
        if(length(res) > 0){
            res <- res[order(padjust)]
        }
    }
    return(res)
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
#' @param subsets A named list of named lists specifying any number of gene 
#'              subsets (can differ per sample). For each subset, FDR correction
#'              will be limited to genes in the subset, and aberrant 
#'              events passing the FDR cutoff will be reported for each subset 
#'              separately. See the examples for how to use this argument. 
#' @param returnTranscriptomewideResults If subsets of genes to test are 
#'              provided (i.e. \code{!is.null(subsets)}), this parameter 
#'              indicates whether additionally the transcritome-wide results 
#'              should be returned as well (default).
#' @param ... Additional arguments, currently not used
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
#' # example of retrieving results with FDR correction limited to a 
#' # set of genes of interest
#' genesOfInterest <- list("sample_1"=sample(rownames(ods), 3), 
#'                          "sample_2"=sample(rownames(ods), 8), 
#'                          "sample_6"=sample(rownames(ods), 5))
#' res_genesOfInterest <- results(ods, 
#'                          subsets=list("Example_subset"=genesOfInterest),
#'                          all=TRUE,
#'                          returnTranscriptomewideResults=FALSE)
#' res_genesOfInterest
#' 
#' @name results
#' @rdname results
#' @aliases results results,OutriderDataSet-method
#' 
#' @export
setMethod("results", "OutriderDataSet", compileResultsAll.OUTRIDER)

