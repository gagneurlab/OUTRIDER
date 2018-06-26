#'
#' OutriderDataSet class and constructors 
#' 
#' @description The OutriderDataSet class is designed to store the whole OUTRIDER data set
#' needed for an analysis. It is a subclass of 
#' \code{RangedSummarizedExperiment}. All calculated values and results are 
#' stored as assays or as annotation in the mcols structure provided by the
#' \code{RangedSummarizedExperiment} class.
#' 
#' @param se An RangedSummarizedExperiment object or any object which inherits
#'             from it and contains a count matrix as the first element in the
#'             assay list.
#' @param countData A simple count matrix. If dim names are provided, they have
#'             to be unique. This is only used if no \code{se} object is
#'             provided.
#' @param colData Additional to the count data a \code{DataFrame} can be 
#'             provided to annotate the samples. 
#' @param ... Further arguments can be passed to
#'             \code{\link[DESeq2]{DESeqDataSet}}, which is used to parse the
#'             user input and create the initial 
#'             \code{RangedSummarizedExperiment} object. 
#' @return An OutriderDataSet object.
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}, 
#'             Felix Brechtmann \email{brechtma@@in.tum.de}
#' 
#' @rdname OutriderDataSet-class
#' 
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet()
#' ods
#'    
setClass("OutriderDataSet", contains="RangedSummarizedExperiment")

#
# check sample annotation within the colData slot of the SE object
#
validateCounts <- function(object) {
    if(!"counts" %in% assayNames(object)){
        return("No counts are detected. Please provide a count matrix.")
    }
    if(!is.integer(assay(object, "counts"))){
        return("Please provide an integer count table.")
    }
    NULL
}

checkNames <- function(object){
    n <- rownames(object)
    if(!is.null(n) && any(duplicated(n))){
        return("Please provide unique rownames or no rownames at all.")
    }
    n <- colnames(object)
    if(!is.null(n) && any(duplicated(n))){
        return("Please provide unique colnames or no colnames at all.")
    }
}


## general validate function
validateOutriderDataSet <- function(object) {
    c(
        checkNames(object),
        validateCounts(object)
    )
}
setValidity("OutriderDataSet", validateOutriderDataSet)


## show method for OutriderDataSet
showOutriderDataSet <- function(object) {
    cat("class: OutriderDataSet\n")
    show(as(object, "SummarizedExperiment"))
}

setMethod("show", "OutriderDataSet", function(object) {
    showOutriderDataSet(object)
})


#' @rdname OutriderDataSet-class
#' @export
OutriderDataSet <- function(se, countData, colData, ...) {
    
    # use SummarizedExperiment object
    if(!missing(se)){
        if(!is(se, "SummarizedExperiment")){
            stop("'se' must be a RangedSummarizedExperiment object")
        }
        se <- DESeqDataSet(se, design=~1, ...)
    
    # use raw count data
    } else if(!missing(countData)){
        if(missing(colData)){
            cols <- colnames(countData)
            if(is.null(cols)){
                cols <- paste0("sample", seq_len(ncol(countData)))
            }
            colData <- DataFrame(sampleID=cols)
            rownames(colData) <- colData[['sampleID']]
            colnames(countData) <- colData[['sampleID']]
        }
        se <- DESeqDataSetFromMatrix(countData=countData, 
                colData=colData, design=~1, ...)
        
    # nothing provided
    } else {
        stop("At least one of the se or countData argument has to be defined!")
    }
    
    obj <- new("OutriderDataSet", se)
    validObject(obj)
    
    return(obj)
}

#'
#' Create example data sets for OUTRIDER
#' 
#' Creates an example data set from a file or simulates a data set based 
#' on random counts following a negative binomial distribution with injected
#' outliers with a fixed z score away from the mean of the gene.
#' 
#' @param dataset If "none", the default, an example data set is simulated. 
#'             One can also use example data set included in the package by
#'             specifying 'GTExSkinSmall' or 'KremerNBaderSmall'
#' @param n Number of simulated genes 
#' @param m Number of simulated samples
#' @param freq Frequency of in-silico outliers
#' @param zScore Absolute z score of in-silico outliers (default 6).
#' @param inj Determines whether counts are injected with the strategy 
#'            ('both', 'low', 'high'), default is 'both'.
#' @param ... Further arguments to \code{makeExampleDESeqDataSet}
#'
#' @return An OutriderDataSet containing an example dataset. Depending on the
#'             parameters it is based on a real data set or it is simulated
#' 
#' @examples
#' # A generic dataset 
#' ods1 <- makeExampleOutriderDataSet()
#' ods1
#' 
#' # A generic dataset with specificed sample size and injection method
#' ods2 <- makeExampleOutriderDataSet(n=500, m=50, inj='low')
#' ods2
#' 
#' # A subset of a real world dataset from GTEx 
#' ods3 <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods3
#' 
#' @export 
makeExampleOutriderDataSet <- function(n=1000, m=100, freq=1E-2, zScore=6, 
                    inj=c('both', 'low', 'high'),
                    dataset=c('none', 'GTExSkinSmall', 'KremerNBaderSmall'),
                    ...){
    
    dataset <- match.arg(dataset)
    inj <- match.arg(inj)
    if(dataset != 'none'){
        file <- system.file("extdata", paste0(dataset, ".tsv"), 
                package="OUTRIDER", mustWork=TRUE)
        countData <- read.table(file)
        return(OutriderDataSet(countData=countData))
    }
    
    ans <- OutriderDataSet(se=makeExampleDESeqDataSet(n=n, m=m*2, ...))
    ans <- ans[,seq_len(m)]
    colnames(ans)[m] <- paste0("sample", m)
    ans <- injectOutliers(ans, freq = freq, zScore = zScore, inj=inj)
    return(ans)
}

