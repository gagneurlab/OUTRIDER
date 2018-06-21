#######################
## OutriderDataSet
## ====================


#' OutriderDataSet
#'
#' This class is designed to store the whole OUTRIDER data set
#' needed for an analysis of a disease cohort
#'
#' @author Christian Mertes \email{mertes@@in.tum.de}
setClass("OutriderDataSet", contains="RangedSummarizedExperiment")

## Validity
## ========

#
# check sample annotation within the colData slot of the SE object
#
validateCounts <- function(object) {
    cts <- assays(object)[['counts']]
    if(!is.null(cts)){
        if(!is.integer(cts)){
            return("Please provide a integer count table.")
        }
    }
    NULL
}

checkRowNames <- function(object){
    NULL
}


## general validate function
validateOutriderDataSet <- function(object) {
    c(
        checkRowNames(object),
        validateCounts(object)
    )
}
setValidity("OutriderDataSet", validateOutriderDataSet)


## Cosmetics (the show function)
## =============================

## show method for OutriderDataSet
showOutriderDataSet <- function(object) {
    cat("class: OutriderDataSet\n")
    show(as(object, "SummarizedExperiment"))
}

setMethod("show", "OutriderDataSet", function(object) {
    showOutriderDataSet(object)
})

## Constructor
## ==========

#'
#' The constructor function for OutriderSettings
#' 
#' Eather a RangedSummarizedExperiment object or a count table has to be 
#' provided.
#' 
#' @param se RangedSummarizedExperiment object
#' @param countData countData 
#' @param colData additional Annotation Data
#' @param ... Any parameters corresponding to the slots and their possible
#' values. See \linkS4class{OutriderDataSet}
#' @return A OutriderDataSet object.
#' @author Christian Mertes \email{mertes@@in.tum.de}
#' 
#' @examples
#'     ods <- makeExampleOutriderDataSet()
#'     ods
#'     
#' @export
OutriderDataSet <- function(se=NULL, countData=NULL, colData=NULL, ...) {
    createdCondition <- FALSE
    if(!is.null(se)){
        if(!is(se, "SummarizedExperiment")){
            stop("'se' must be a RangedSummarizedExperiment object")
        }
        if(!"condition" %in% colnames(colData(se))){
            createdCondition <- TRUE
            colData(se)$condition <- factor(paste0('C', 1:dim(se)[2]))
        }
        se <- DESeqDataSet(se, design=~condition, ...)
        
    } else if(!is.null(countData)){
        if(is.null(colData)){
            colData <- DataFrame(
                condition = factor(paste0("C", 1:dim(countData)[2])))
            rownames(colData) <- paste0("sample", 1:dim(countData)[2])
            if(!is.null(colnames(countData))){
                rownames(colData) <- colnames(countData)
            }
        } else {
            if(!"condition" %in% colnames(colData)){
                createdCondition <- TRUE
                colData$condition <- factor(paste0("C", 1:dim(countData[2])))
            }
        }
        se <- DESeqDataSetFromMatrix(countData=countData, 
                colData=colData, design=~condition, ...)
        
    } else {
        stop("One of the se or countData parameter has to be defined!")
    }
    
    obj <- new("OutriderDataSet", se)
    validObject(obj)
    
    if(createdCondition){
        colData(obj)$condition <- NULL
    }
    
    return(obj)
}

#'
#' Create example data sets for OUTRIDER
#' 
#' Creates an example data set from a file or a generates a random counts.
#' 
#' @param dataset here one can select from the two example data sets.
#'             One of 'none', 'GTExSkinSmall', or 'KremerNBaderSmall'.
#' @param n number of simulated genes 
#' @param m number of simulated samples.
#' @param ... further arguments to \code{makeExampleDESeqDataSet}
#' @return An OutriderDataSet containing an example dataset. Depending on the
#'             parameters it is based on a real data set or on simulated counts.
#' 
#' @examples
#' # A generic dataset 
#' ods1 <- makeExampleOutriderDataSet()
#' ods1
#' 
#' # A generic dataset with specificed sample size
#' ods2 <- makeExampleOutriderDataSet(n=10, m=100)
#' ods2
#' 
#' # A subset of a real world dataset from GTEx 
#' ods3 <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods3
#' 
#' @export 
makeExampleOutriderDataSet <- function(n=1000, m=100, freq = 1E-2, zScore = 6, inj='both', ..., 
                    dataset=c('none', 'GTExSkinSmall', 'KremerNBaderSmall')){
    dataset <- match.arg(dataset)
    if(dataset != 'none'){
        file <- system.file("extdata", paste0(dataset, ".tsv"), 
                package="OUTRIDER", mustWork=TRUE)
        countData <- read.table(file)
        return(OutriderDataSet(countData=countData))
    }
    
    ans <- OutriderDataSet(se=makeExampleDESeqDataSet(n=n, m=m*2, ...))
    ans <- ans[,c(1:m-1, m+1)]
    colnames(ans)[m] <- paste0("sample", m)
    ans <- injectOutliers(ans, freq = freq, zScore = zScore, inj=inj)
    return(ans)
}

