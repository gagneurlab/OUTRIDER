filterExpression.OUTRIDER <- function(object, gtfFile, fpkmCutoff=1, 
                    percentile=0.95, filterGenes=TRUE, savefpkm=FALSE, 
                    minCounts=FALSE, addExpressedGenes=TRUE, ...){
    object <- filterMinCounts(object, filterGenes=filterGenes, 
            verbose=minCounts)
    if(isTRUE(minCounts)){
        return(object)
    }
    if(!missing(gtfFile)){
        object <- computeGeneLength(object, gtfFile=gtfFile, ...)
    }
    filterExp(object, fpkmCutoff=fpkmCutoff, percentile=percentile,
            filterGenes=filterGenes, savefpkm=savefpkm, 
            addExpressedGenes=addExpressedGenes)
}

#' 
#' Filter expression
#' 
#' To filter out non expressed genes this method uses the FPKM values to 
#' get a comparable value over genes. For each gene, if the pth-
#' \code{percentile} is greater than the \code{fpkmCutoff} value, it passes the
#' filter. To calcute the FPKM values the user needs to provide a GTF file or 
#' the basepair parameter as described in \code{\link[DESeq2]{fpkm}}.
#' 
#' @rdname filterExpression
#' @param object An OutriderDataSet object
#' @param filterGenes If TRUE, the default, the object is subseted.
#' @param minCounts If TRUE, only genes with 0 counts in all samples are 
#'             filtered
#' @param gtfFile A txDb object or a GTF/GFF file to be used as annotation
#' @param fpkmCutoff The threshold for filtering based on the FPKM value
#' @param percentile a numeric indicating the percentile FPKM value to compare
#'             against the \code{fpkmCutoff}
#' @param savefpkm If TRUE, the FPKM values are saved as assay
#' @param addExpressedGenes If TRUE (default), adds 5 columns to the
#'             \code{colData} with information regarding the number of 
#'             expressed genes per sample
#' @param ... Additional arguments passed to \code{computeGeneLength}
#' @return An OutriderDataSet containing the \code{passedFilter} column, which
#'             indicates if the given gene passed the filtering threshold. If
#'             \code{filterGenes} is TRUE the object is already subsetted.
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' annotationFile <- system.file("extdata", 
#'     "gencode.v19.genes.small.gtf.gz", package="OUTRIDER")
#' ods <- filterExpression(ods, annotationFile, filterGenes=FALSE)
#' 
#' mcols(ods)['passedFilter']
#' fpkm(ods)[1:10,1:10]
#' dim(ods)
#' 
#' ods <- ods[mcols(ods)[['passedFilter']]]
#' dim(ods)
#' 
#' @export
setMethod("filterExpression", "OutriderDataSet", filterExpression.OUTRIDER)

filterExp <- function(ods, fpkmCutoff, percentile, filterGenes, savefpkm, 
                    addExpressedGenes){
    fpkm <- fpkm(ods)
    if(savefpkm){
        assay(ods, 'fpkm', withDimnames=FALSE) <- fpkm
    }
    passed <- rowQuantiles(fpkm, probs=percentile) > fpkmCutoff
    mcols(ods)['passedFilter'] <- passed
    
    if(addExpressedGenes == TRUE){
        dt <- computeExpressedGenes(fpkm, cutoff=fpkmCutoff, 
                percentile=percentile)
        goodCols <- !colnames(colData(ods)) %in% colnames(dt[,-1])
        colData(ods) <- DataFrame(row.names=colData(ods)$sampleID,
                merge(colData(ods)[,goodCols, drop=FALSE], 
                        dt, by="sampleID", sort=FALSE))
    }
    
    message(paste0(sum(!passed), ifelse(filterGenes,
            " genes are filtered out. ", " genes did not pass the filter. "),
            "This is ", signif(sum(!passed)/length(passed)*100, 3), 
            "% of the genes."))
    if(filterGenes==TRUE){
        ods <- ods[passed == TRUE,]
    }
    
    validObject(ods)
    return(ods)
}

#' 
#' Extracting the gene length from annotations
#' 
#' Computes the length for each gene based on the given GTF file or annotation.
#' Here the length of a gene is defind by the total number of bases covered
#' by exons.
#' 
#' @param ods An OutriderDataSet for which the gene length should be computed.
#' @param gtfFile Can be a GTF file or an txDb object with annotation.
#' @param format The format parameter from \code{makeTxDbFromGFF}
#' @param mapping If set, it is used to map gene names between the GFT and the
#'             ods object. This should be a 2 column data.frame: 
#'             1. column GTF names and 2. column ods names.
#' @param ... further arguments to \code{makeTxDbFromGFF}
#' 
#' @return An OutriderDataSet containing a \code{basepairs} column with the 
#'             calculated gene length. Accessable through 
#'             \code{mcols(ods)['baisepairs']}
#' 
#' @examples 
#' 
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' annotationFile <- system.file("extdata", "gencode.v19.genes.small.gtf.gz",
#'         package="OUTRIDER")
#' ods <- computeGeneLength(ods, annotationFile)
#' 
#' mcols(ods)['basepairs']
#' fpkm(ods)[1:10,1:10]
#' 
#' @export
computeGeneLength <- function(ods, gtfFile, format='gtf', mapping=NULL, ...){
    checkOutriderDataSet(ods)
    
    if(!is(gtfFile, "TxDb")){
        #created txdb object
        txdb <- makeTxDbFromGFF(gtfFile, format=format, ...)
    } else {
        txdb <- gtfFile
    }
    
    # collect the exons per gene id
    exons_list_per_gene <- exonsBy(txdb,by="gene")
    
    # for each gene, reduce all the exons to a set of non overlapping exons,
    # calculate their lengths (widths) and sum then
    widths <- width(reduce(exons_list_per_gene))
    totalexonlength <- vapply(widths, sum, numeric(1))
    if(!is.null(mapping)){
        if(is.data.table(mapping)){
            mapping <- as.data.frame(mapping)
        }
        if(is.data.frame(mapping)){
            tmpmapping <- mapping[,2]
            names(tmpmapping) <- mapping[,1]
            mapping <- tmpmapping
        }
        names(totalexonlength) <- mapping[names(totalexonlength)]
    }
    mcols(ods)['basepairs'] <- totalexonlength[rownames(ods)]
    
    # checking for NAs in basepair annotation
    if(any(is.na(mcols(ods)['basepairs']))){
        missingNames <- rownames(ods)[is.na(mcols(ods)['basepairs'])]
        warning(paste0("Some genes (n=", length(missingNames), 
                ") are not found in the annotation. Setting 'basepairs' == 1. ",
                "The first 6 names are:\n", paste(collapse=", ",
                        missingNames[seq_len(min(6, length(missingNames)))])))
        mcols(ods[is.na(mcols(ods)['basepairs'])[,1]])['basepairs'] <- 1
    }
    
    validObject(ods)
    return(ods)
}


#' 
#' Calculate FPM and FPKM values
#' 
#' This is the fpm and fpkm function from DESeq2. For more details see: 
#' \code{\link[DESeq2]{fpkm}} and \code{\link[DESeq2]{fpm}}
#' 
#' @name fpkm
#' @rdname fpkm
#' @aliases fpkm fpm
#' @seealso \code{\link[DESeq2]{fpkm}} \code{\link[DESeq2]{fpm}}
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' mcols(ods)['basepairs'] <- round(rnorm(nrow(ods), 1000, 500))
#' 
#' mcols(ods)['basepairs']
#' fpkm(ods)[1:10,1:10]
#' fpm(ods)[1:10,1:10]
#' 
#' @export fpkm
#' @export fpm
NULL

#'
#' Filter zero counts
#' 
#' @noRd
filterMinCounts <- function(x, filterGenes=FALSE, verbose=TRUE){
    passed <- !checkCountRequirements(x, test=TRUE)
    mcols(x)['passedFilter'] <- passed
    
    if(isTRUE(verbose)){
        message(paste0(sum(!passed), " genes did not pass the filter due to ", 
                "zero counts. This is ", 
                signif(sum(!passed)/length(passed)*100, 3), "% of the genes."))
    }
    
    if(isTRUE(filterGenes)){
        x <- x[passed,]
    }
    return(x)
}


#'
#' Returns a data.table with summary statistics per sample on number of genes
#' 
#' @noRd
computeExpressedGenes <- function(fpkm, cutoff=1, percentile=0.95){
    
    # Get cells that passed cutoff
    cutoffPassedMatrix <- fpkm > cutoff
    
    # Remove rows where no genes passed the cutoff
    cutoffPassedMatrix <- cutoffPassedMatrix[rowSums(cutoffPassedMatrix) > 0,]
    
    # Make a data.table with the expressed genes
    expGenesDt <- data.table(sampleID = colnames(cutoffPassedMatrix), 
            expressedGenes = colSums(cutoffPassedMatrix))
    
    # order by sample rank
    setorder(expGenesDt, "expressedGenes")
    cutoffPassedMatrix <- cutoffPassedMatrix[, expGenesDt$sample]
    
    # Get the cummulative sum of expressed genes
    cumSumMatrix <- rowCumsums(cutoffPassedMatrix + 0)
    
    # Get the genes that appear in at least 1 sample
    expGenesDt$unionExpressedGenes <- colSums(cumSumMatrix > 0)
    
    # Get the genes common in all samples
    expGenesDt$intersectionExpressedGenes <- rowSums(
            t(cumSumMatrix) == seq_col(cumSumMatrix))
    
    # Get number of genes passing the filtering
    expGenesDt$passedFilterGenes <- colSums(
            t(t(cumSumMatrix)/seq_col(cumSumMatrix)) >= 1-percentile)
    
    # Rank for plotting
    expGenesDt[, c("expressedGenesRank"):=list(.I)]
    
    return(expGenesDt)
}
