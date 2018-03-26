#' 
#' Filter expression
#' 
#' To filter out none expressed genes this method uses the fpkm values to 
#' get a comparable value over genes. To calcute the fpkm values the user 
#' needs to provide a GTF file or the basepair parameter as described in 
#' \code{\link[DESeq2]{fpkm}}.
#' 
#' @rdname filterExpression
#' @param x An OutriderDataSet object
#' @param filterGenes if TRUE, the default, the object is subseted.
#' @param onlyZeros filter only based on zero counts on a gene
#' @param gtfFile a txdb object or a GTF/GFF file to be used as annotation
#' @param fpkmCutoff the threshold for filtering based on the FPKM value
#' @param savefpkm if TRUE the FPKM values are saved as assay
#' @param ... additional arguments passed to \code{computeGeneLength}
#' @return An OutriderDataSet containing the \code{passedFilter} column, which
#'             indicates if the given gene passed the filtering threshold. If
#'             filterGenes is TRUE the object is already subsetted.
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' annotationFile <- system.file("extdata", 
#'     "gencode.v19.genes.small.gtf.gz", package="OUTRIDER")
#' ods <- filterExpression(ods, annotationFile)
#' 
#' mcols(ods)['passedFilter']
#' fpkm(ods)[1:10,1:10]
#' dim(ods)
#' 
#' ods <- filterExpression(ods, annotationFile, filterGenes=TRUE)
#' dim(ods)
#' 
#' @export
setGeneric("filterExpression", 
        function(x, ...) standardGeneric("filterExpression"))

#' @rdname filterExpression
#' @export
setMethod("filterExpression", "OutriderDataSet", function(x, gtfFile=NULL, 
                    fpkmCutoff=1, filterGenes=FALSE, savefpkm=FALSE, 
                    onlyZeros=FALSE, ...){
    x <- filterZeros(x, filterGenes=filterGenes)
    if(onlyZeros == TRUE){
        return(x)
    }
    if(!is.null(gtfFile)){
        x <- computeGeneLength(x, gtfFile=gtfFile, ...)
    }
    filterExp(x, fpkmCutoff=fpkmCutoff, filterGenes=filterGenes, 
            savefpkm=savefpkm)
})

filterExp <- function(ods, fpkmCutoff=1, filterGenes=filterGenes, 
                savefpkm=savefpkm){
    fpkm <- fpkm(ods)
    if(savefpkm){
        assays(ods)[['fpkm']]<-fpkm
    }
    passed <- rowQuantiles(fpkm, probs=c(0.95)) > fpkmCutoff

    mcols(ods)['passedFilter'] <- passed
    validObject(ods)
    
    if(filterGenes==TRUE){
        ods <- ods[passed == TRUE]
    }
    message(paste0(sum(passed), ifelse(filterGenes, 
            " genes are filtered out. ", " genes did not passed the filter. "), 
            "This is ", signif(sum(passed)/length(passed)*100, 3), 
            "% of the genes."))
    return(ods)
}

#' 
#' Computes for each gene based on the GTF file the exon length
#' 
#' @param ods An OutriderDataSet for which the gene length should be computed.
#' @param gtfFile Can be a gft file or an txDb object with annotation.
#' @param format the format parameter from \code{makeTxDbFromGFF}
#' @param mapping if set, it is used to map gene names between gtf and ods. 
#'             This should be a 2 column data.frame: 1. column GTF names 
#'             and 2. column ods names.  
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
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' ods <- makeExampleOutriderDataSet(dataset="KremerNBaderSmall")
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' mapping <- select(org.Hs.eg.db, keys=keys(txdb, keytype = "GENEID"), 
#'         keytype="ENTREZID", columns=c("SYMBOL"))
#' ods <- computeGeneLength(ods, txdb, mapping=mapping)
#' 
#' mcols(ods)['basepairs']
#' fpkm(ods)[1:10,1:10]
#' 
#' @export
computeGeneLength <- function(ods, gtfFile, format='gtf', mapping=NULL, ...){
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
                "The first 6 names are:\n", 
                paste(missingNames[1:min(6, length(missingNames))], 
                        collapse=", ")))
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
#' @inheritParams DESeq2::fpm
#' @inheritParams DESeq2::fpkm
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' annotationFile <- system.file("extdata", 
#'     "gencode.v19.genes.small.gtf.gz", package="OUTRIDER")
#' ods <- computeGeneLength(ods, annotationFile)
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
filterZeros <- function(x, filterGenes=FALSE){
    stopifnot(is(x, 'OutriderDataSet'))
    passed <- rowMax(counts(x)) > 0
    mcols(x)['passedFilter'] <- passed
    
    message(paste0(sum(passed), " genes ", ifelse(filterGenes, 
            "are filtered out ", "did not passed the filter "), 
            "due to zero counts. This is ", 
            signif(sum(passed)/length(passed)*100, 3), 
            "% of the genes."))
    
    if(isTRUE(filterGenes)){
        x <- x[passed,]
    }
    return(x)
}
