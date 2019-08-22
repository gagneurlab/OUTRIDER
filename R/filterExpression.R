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
#' @param x An OutriderDataSet object
#' @param filterGenes If TRUE, the default, the object is subseted.
#' @param minCounts If TRUE, only genes with 0 counts in all samples are 
#'             filtered
#' @param gtfFile A txDb object or a GTF/GFF file to be used as annotation
#' @param fpkmCutoff The threshold for filtering based on the FPKM value
#' @param percentile a numeric indicating the percentile FPKM value to compare
#'             against the \code{fpkmCutoff}
#' @param savefpkm If TRUE, the FPKM values are saved as assay
#' @param addExpressedGenes If TRUE, adds 5 columns to the \code{colData} with
#'            information regarding the number of expressed genes per sample
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
setGeneric("filterExpression", 
           function(x, ...) standardGeneric("filterExpression"))

filterExpression.OUTRIDER <- function(x, gtfFile, fpkmCutoff=1, 
                                      percentile=0.95, filterGenes=TRUE, savefpkm=FALSE, 
                                      minCounts=FALSE, addExpressedGenes=FALSE, ...){
  x <- filterMinCounts(x, filterGenes=filterGenes)
  if(minCounts == TRUE){
    return(x)
  }
  if(!missing(gtfFile)){
    x <- computeGeneLength(x, gtfFile=gtfFile, ...)
  }
  filterExp(x, fpkmCutoff=fpkmCutoff, percentile=percentile,
            filterGenes=filterGenes, savefpkm=savefpkm, 
            addExpressedGenes=addExpressedGenes)
}

#' @rdname filterExpression
#' @export
setMethod("filterExpression", "OutriderDataSet", filterExpression.OUTRIDER)

filterExp <- function(ods, fpkmCutoff=1, percentile=0.95,
                      filterGenes=filterGenes, savefpkm=savefpkm, 
                      addExpressedGenes=FALSE){
  fpkm <- fpkm(ods)
  if(savefpkm){
    assays(ods)[['fpkm']]<-fpkm
  }
  passed <- rowQuantiles(fpkm, probs=percentile) > fpkmCutoff
  
  mcols(ods)['passedFilter'] <- passed
  
  if(addExpressedGenes == TRUE){
    dt <- computeExpressedGenes(fpkm, cutoff = fpkmCutoff, percentile = percentile)
    colData(ods) <- merge(colData(ods), dt, sort = FALSE)
  }
  validObject(ods)
  
  if(filterGenes==TRUE){
    ods <- ods[passed == TRUE]
  }
  message(paste0(sum(!passed), ifelse(filterGenes, 
                                      " genes are filtered out. ", " genes did not passed the filter. "), 
                 "This is ", signif(sum(!passed)/length(passed)*100, 3), 
                 "% of the genes."))
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
                   "The first 6 names are:\n", 
                   paste(missingNames[seq_len(min(6, length(missingNames)))],
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
filterMinCounts <- function(x, filterGenes=FALSE){
  passed <- !checkCountRequirements(x, test=TRUE)
  mcols(x)['passedFilter'] <- passed
  
  message(paste0(sum(!passed), " genes did not passed the filter due to ", 
                 "zero counts. This is ", signif(sum(!passed)/length(passed)*100, 3),
                 "% of the genes."))
  
  if(isTRUE(filterGenes)){
    x <- x[passed,]
  }
  return(x)
}


#'
#' Returns a data.table with the N of expressed genes per sample
#' 
#' @noRd
computeExpressedGenes <- function(x, cutoff=1, percentile=0.95){
  
  # Get cells that passed cutoff
  logic_mat <- x > cutoff
  
  # Remove rows where no genes passed the cutoff
  logic_mat <- logic_mat[rowSums(logic_mat) > 0, ]
  
  # Make a data.table with the expressed genes
  dt <- data.table(sampleID = colnames(logic_mat), expressedGenes = colSums(logic_mat))
  setorder(dt, expressedGenes)
  
  logic_mat <- logic_mat[, dt$sample]
  
  # Get the cummulative sum of expressed genes
  cum_sum_mat <- t(apply(logic_mat, 1, cumsum))
  
  # Get the genes common in all samples
  dt$unionExpressedGenes <- colSums(cum_sum_mat > 0)
  
  # Get the genes that appear in at least 1 sample
  dt$intersectionExpressedGenes <- colSums(sapply(1:ncol(cum_sum_mat), function(j) cum_sum_mat[,j] == j))
  
  dt$passedFilterGenes <- colSums(t(t(cum_sum_mat) / 1:ncol(cum_sum_mat)) >= 1-percentile)
  
  # Rank for plotting
  dt[, expressedGenesRank := .I]
  
  return(dt)
}
