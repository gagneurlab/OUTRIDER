compileResults.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                    round=2, all=FALSE, subsetName=NULL, ...){
    
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
    
    meanCorrectedCounts <- rowMeans(counts(object, normalized=TRUE))
    meanRawCounts       <- rowMeans(counts(object, normalized=FALSE))
    if(isFALSE(all)){
        abByGene <- aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="gene",
                subsetName=subsetName)
        abBySample <- aberrant(object, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, by="sample",
                subsetName=subsetName)
        object <- object[abByGene > 0, abBySample > 0]
        meanCorrectedCounts <- meanCorrectedCounts[abByGene > 0]
        meanRawCounts <- meanRawCounts[abByGene > 0]
    }
    
    # warning if no rows to return 
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
            meanRawcounts=NA_real_,
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
        meanRawcounts    = meanRawCounts,
        normcounts       = c(counts(object, normalized=TRUE)),
        meanCorrected    = meanCorrectedCounts,
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
        padj_rank        = c(apply(padj(object, subsetName=subsetName), 2, rank)),
        FDR_set          = ifelse(is.null(subsetName), "transcriptome-wide", 
                                    subsetName)
    )
    
    # round columns if requested
    if(is.numeric(round)){
        devNull <- lapply(c("normcounts", "zScore", "l2fc", "theta", 
                "meanRawcounts", "meanCorrected"),
                function(x) ans[,c(x):=round(get(x), as.integer(round))] )
    }
    
    # 
    # keep only aberrant events and sort by padj value
    # 
    if(isFALSE(all)){
        ans <- ans[aberrant == TRUE]
        
    } else{
        # if return full subset is requested, retrieve those as all non-NA padj
        ans <- ans[!is.na(padjust),]
    }
    ans <- ans[order(padjust)]
    
    return(ans)
}

compileResultsAll.OUTRIDER <- function(object, padjCutoff=0.05, zScoreCutoff=0, 
                                    round=2, all=FALSE, 
                                    returnTranscriptomewideResults=TRUE, ...){
    #
    # input check and parsing
    # 
    checkOutriderDataSet(object)
    checkFullAnalysis(object)
    
    # get all padjust assays (transcriptome-wide and on subsets)
    padjAssays <- grep('padjust', assayNames(object), value=TRUE)
    
    # dont retrieve transcriptome-wide results if requested
    if(isFALSE(returnTranscriptomewideResults)){
        if(length(padjAssays) == 1 && padjAssays == 'padjust'){
            warning("Retrieving transcriptome-wide results as no other ",
                    "padjust assays are available in the ods object.")
        } else{
            padjAssays <- padjAssays[padjAssays != 'padjust']
        }
    }
    
    # get results for all available padjust assays
    resSubsets <- lapply(padjAssays, function(padjAssay){
        if(padjAssay == 'padjust'){
            subsetName <- NULL
        } else{
            subsetName <- gsub("padjust_", "", padjAssay)
        }
        resSub <- compileResults.OUTRIDER(object=object, 
                            padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, 
                            all=all, round=round, 
                            subsetName=subsetName, ...)
        return(resSub)
    })
    res <- rbindlist(resSubsets)
    
    # sort it if existing
    if(nrow(res) > 0){
        # get aberrant columns from transcriptome-wide results
        res_tw <- res[FDR_set == "transcriptome-wide", 
                        .(geneID, sampleID, padj_rank, aberrant, 
                            AberrantBySample, AberrantByGene)] 
        
        # dcast to have only one row per gene-sample combination
        res <- dcast(res[, !c("padj_rank", "aberrant", "AberrantBySample", 
                            "AberrantByGene"), with=FALSE], 
                    ... ~ FDR_set, value.var="padjust")
        
        # merge again with aberrant columns for transcriptome-wide results
        if(isTRUE(returnTranscriptomewideResults)){
            # rename column back to padjust for tw results after dcast
            if("transcriptome-wide" %in% colnames(res)){
                setnames(res, "transcriptome-wide", "padjust")
            } else{
                res[, padjust := NA]
            }
            
            # merge with aberrant columns
            res <- merge(res, res_tw, by=c("geneID", "sampleID"), all.x=TRUE)
            
            # set aberrant to FALSE if only significant in subset results
            res[is.na(aberrant), aberrant := FALSE]
            res[is.na(padjust), padjust := padj(object)[cbind(geneID,sampleID)]]
            
            # restore previous column order
            setcolorder(res, colnames(resSubsets[[1]][,!"FDR_set", with=FALSE]))
            
            # order rows based on transcriptome-wide FDR
            res <- res[order(padjust)]
        }
        
        # rename padjust columns for other gene sets to padjust_setname
        for(subsetAssayName in padjAssays[padjAssays != "padjust"]){
            setnames(res, gsub("padjust_", "", subsetAssayName), 
                        subsetAssayName, skip_absent=TRUE)
        }
    } else{
        res <- res[, !"FDR_set", with=FALSE]
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
#'             data.table of length samples x genes. 
#' @param returnTranscriptomewideResults If FDR corrected pvalues for subsets 
#'              of genes of interest have been calculated, this parameter 
#'              indicates whether additionally the transcriptome-wide results 
#'              should be returned as well (default), or whether only results 
#'              for those subsets should be retrieved.
#' @param ... Additional arguments, currently not used
#' @return A data.table where each row is an outlier event and the columns
#'    contain additional information about this event. In details the 
#'    table contains: 
#'    \item{sampleID/geneID}{The gene or sample ID as provided by the 
#'          user, e.g. \code{rowData(ods)} and \code{colData(ods)},
#'          respectively.}
#'    \item{pValue/padjust}{The nominal P-value and the FDR corrected
#'          P-value (transcriptome-wide) indicating the outlier status.}
#'    \item{zScore/l2fc}{The z score and log\eqn{_2}{[2]} fold change 
#'        as computed by \code{\link[OUTRIDER]{computeZscores}}.}
#'    \item{rawcounts}{The observed read counts.}
#'    \item{normcounts}{The expected count given the fitted 
#'          autoencoder model for the given gene-sample combination.}
#'    \item{meanRawcounts/meanCorrected}{For this gene, the mean of the 
#'         observed or expected counts, respectively, given the fitted 
#'         autoencoder model.}
#'    \item{theta}{The dispersion parameter of the NB distribution 
#'          for the given gene.}
#'    \item{aberrant}{The transcriptome-wide outlier status of this event: 
#'          \code{TRUE} or \code{FALSE}.}
#'    \item{AberrantBySample/AberrantByGene}{Number of outliers for the 
#'          given sample or gene (transcriptome-wide), respectively.}
#'    \item{padj_rank}{Rank of this outlier event within the given sample.}
#'    \item{padjust_FDRset}{The FDR corrected P-value with respect to the 
#'          gene subset called 'FDRset', if gene subsets were specified 
#'          during the P-value computation. Find more details at 
#'          \code{\link[OUTRIDER]{computePvalues}}.}
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
#' genesOfInterest
#' ods <- computePvalues(ods, subsets=list("exampleSubset"=genesOfInterest))
#' res <- results(ods, all=TRUE, returnTranscriptomewideResults=FALSE)
#' res
#' 
#' @name results
#' @rdname results
#' @aliases results results,OutriderDataSet-method
#' 
#' @export
setMethod("results", "OutriderDataSet", compileResultsAll.OUTRIDER)

