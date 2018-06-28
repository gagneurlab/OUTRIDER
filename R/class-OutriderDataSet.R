#'
#' OutriderDataSet class and constructors 
#' 
#' @description The OutriderDataSet class is designed to store the whole 
#' OUTRIDER data set needed for an analysis. It is a subclass of 
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
    show(as(object, "RangedSummarizedExperiment"))
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
            stop("'se' must be a SummarizedExperiment object")
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
#' @param sizeFactors Artificial Size Factors
#' @param betaSD Standard deviation.
#' @param dispMeanRel Mean dispersion relation ship. 
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
#' ods2 <- makeExampleOutriderDataSet(n=200, m=50, inj='low')
#' ods2
#' 
#' # A subset of a real world dataset from GTEx 
#' ods3 <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods3
#' 
#' @export 
makeExampleOutriderDataSet <- function(n=200, m=80, freq=1E-2, zScore=6, 
                    inj=c('both', 'low', 'high'), sizeFactors = rep(1, m),
                    betaSD=4, dispMeanRel = function(x) 4/x + 0.1,
                    dataset=c('none', 'GTExSkinSmall', 'KremerNBaderSmall')
                    ){
    # load example data set 
    dataset <- match.arg(dataset)
    inj <- match.arg(inj)
    if(dataset != 'none'){
        file <- system.file("extdata", paste0(dataset, ".tsv"), 
                package="OUTRIDER", mustWork=TRUE)
        countData <- read.table(file)
        return(OutriderDataSet(countData=countData))
    }
    
    # generate in-silico data set
    x <- rep(1, m)
    beta <- rnorm(n, 5, betaSD)
    dispersion <- dispMeanRel(2^beta)
    mu <- t(2^(x %*% t(beta)) * sizeFactors)
    countData <- matrix(rnbinom(m * n, mu = mu, size = 1/dispersion), 
            ncol = m)
    
    ## generate in-silico outliers.
    # generate index of injected counts
    index <- matrix(sample(c(0,1,-1), n*m, prob = c(1 - freq, freq/2, freq/2), 
            replace = TRUE), nrow = n)
    # inject on low, high or both sides
    if(inj=='low'){
        index <- -abs(index)
    }
    if(inj=='high'){
        index <- abs(index)
    }
    list_index <- which(index != 0, arr.ind = TRUE)
    
    max_out <- 1E2 * max(countData)
    n_rejected <- 0
    for(i in seq_len(nrow(list_index))){
        row <- list_index[i,'row']
        col <- list_index[i,'col']
        fc <- zScore * sdLogScale(mu[row,col], dispersion[row])
        clcount <- index[row,col]*fc + 
            log2(1 + countData[row,col])
        #multiply size factor again
        art_out <- round(sizeFactors[col]*2^clcount)
        if(art_out < max_out){
            countData[row,col] <- art_out
        }else{
            #remove super large outliers
            index[row,col] <- 0 
            n_rejected <- n_rejected + 1
        }
        
    }
    
    mode(countData) <- "integer"
    
    colnames(countData) <- paste("sample", seq_len(m), sep = "")
    rowRanges <- GRanges("1", IRanges(start=seq_len(n) * 100, width=100))
    names(rowRanges) <- paste0("gene", seq_len(n))
    
    object <- OutriderDataSet(countData = countData, rowRanges = rowRanges)
    trueVals <- DataFrame(trueBeta = beta, trueDisp = dispersion,
            trueMean = rowMeans2(mu))
    mcols(trueVals) <- DataFrame(type = rep("input", ncol(trueVals)), 
            description = c("simulated beta values", "simulated means",
                "simulated dispersion values"))
    mcols(object) <- cbind(mcols(object), trueVals)
    metadata(object)[['trueOutliers']] <- index
    
    return(object)
    
    return(object)
}

#' Aproximates standard deviation of counts in log2 space.
#'
#'@noRd
sdLogScale <- function(mu, disp){
    sqrt(mu*(1+mu/disp))/((mu + 1)*log(2))
}
