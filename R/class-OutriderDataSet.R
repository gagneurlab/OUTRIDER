#' 
#' OutriderDataSet class and constructors 
#' 
#' The OutriderDataSet class is designed to store the whole 
#' OUTRIDER data set needed for an analysis. It is a subclass of 
#' \code{RangedSummarizedExperiment}. All calculated values and results are 
#' stored as assays or as annotation in the mcols structure provided by the
#' \code{RangedSummarizedExperiment} class.
#' 
#' @param se A RangedSummarizedExperiment object or any object which inherits
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
#' ods <- makeExampleOutriderDataSet(dataset="Kremer")
#' ods
#'    
setClass("OutriderDataSet", contains="Outrider2DataSet",
         prototype = list(
             distribution    = "Negative-Binomial",
             preprocessing   = "None",
             transformation  = "sf-log",
             fitModel             = "Autoencoder"
         ))

#' check sample annotation within the colData slot of the SE object
#' @noRd
validateCounts <- function(object) {
    if(!"counts" %in% assayNames(object)){
        return("No counts are detected. Please provide a count matrix.")
    }
    if(!is.integer(assay(object, "counts"))){
        return("Please provide an integer count table.")
    }
    if(min(assay(object, "counts")) < 0){
        return("Please provide a count table with non-negative integers.")
    }
    NULL
}


#' general validate function
#' @noRd
validateOutriderDataSet <- function(object) {
    c(
        checkNames(object),
        validateCounts(object),
        validateModelParameters(object)
    )
}
setValidity("OutriderDataSet", validateOutriderDataSet)


#' show method for OutriderDataSet
#' @noRd
showOutriderDataSet <- function(object) {
    cat("class: OutriderDataSet\n")
    # show(as(object, "RangedSummarizedExperiment"))
    show(as(object, "Outrider2DataSet"))
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
    
    # set sampleID and colnames 
    if(!"sampleID" %in% colnames(se@colData)){
        warning("No sampleID was specified. We will generate a generic one.")
        if(!is.null(colnames(se)) && !any(duplicated(colnames(se)))){
            se@colData[["sampleID"]] <- colnames(se)
        } else {
            se@colData[["sampleID"]] <- paste0("sample_", seq_col(se))
        }
    }
    colnames(se) <- se@colData[["sampleID"]]
    
    obj <- new("OutriderDataSet", se)
    metadata(obj)[["version"]] <- packageVersion("OUTRIDER")
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
#' @param sf Artificial Size Factors
#' @param q number of simulated latend variables.
#' 
#' @return An OutriderDataSet containing an example dataset. Depending on the
#'            parameters it is based on a real data set or it is simulated
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
makeExampleOutriderDataSet <- function(n=200, m=80, q=10, freq=1E-3, zScore=6,
                    inj=c('both', 'low', 'high'), sf=rnorm(m, mean=1, sd=0.1),
                    dataset=c('none', 'GTExSkinSmall', 'KremerNBaderSmall')){
    # load example data set 
    dataset <- match.arg(dataset)
    if(dataset != 'none'){
        file <- system.file("extdata", paste0(dataset, ".tsv"), 
                package="OUTRIDER", mustWork=TRUE)
        countData <- read.table(file)
        return(OutriderDataSet(countData=countData))
    }
    
    theta <- rlnorm(n, meanlog=log(180), sdlog=2)  # dispersion
    logMean <- 5                                   # log offset for mean expres
    sdVec <- rep(0.5, m)                           # sd for H matrix
    
    #
    # Simulate covariates.
    #
    H_true <- matrix(rnorm(m*q), nrow=m, ncol=q)
    D_true <- matrix(rnorm(n*q, sd=sdVec), nrow=n, ncol=q)
    y_true <- D_true %*% t(cbind(H_true))
    mu     <- t(t(exp(rnorm(n, logMean) + y_true))*sf)
    
    # scale it up to overcome the estimation error
    true_sd <- sdLogScale(rowMeans(mu), theta)
    
    #
    # Simulate count Matrix with specified means.
    #
    k <- matrix(rnbinom(m*n, mu=mu, size=theta), nrow=n, ncol=m)
    mode(k) <- 'integer'
    
    #
    # Create Outrider data set
    #
    countData <- DataFrame(k, row.names=paste0("feature_", seq_len(n)))
    colnames(countData) <- paste0("sample_", seq_len(m))
    ods <- OutriderDataSet(countData=countData)
    
    assay(ods, "trueMean", withDimnames=FALSE) <- mu
    assay(ods, "trueSd", withDimnames=FALSE) <- matrix(true_sd, nrow=n, ncol=m)
    mcols(ods)[,"trueTheta"]                 <- theta
    colData(ods)[['trueSizeFactor']]         <- sf
    metadata(ods)[['optimalEncDim']]         <- q
    metadata(ods)[['encDimTable']]           <- data.table(
        encodingDimension=q, evaluationLoss=1, evalMethod='simulation')
    
    #
    # inject outliers
    #
    indexOut <- matrix(nrow=n, ncol=m,
        sample(c(-1,1,0), m*n, replace=TRUE, prob=c(freq/2, freq/2, 1-freq)))
    indexOut <- switch(match.arg(inj),
        low  = -abs(indexOut),
        high =  abs(indexOut),
        indexOut
    )
    normtable <- t(t(k)/sf)
    datasd <- assay(ods, "trueSd")
    lmu <- log2(assay(ods, "trueMean"))
    
    # inject outliers
    max_out <- 1E2 * min(max(k), .Machine$integer.max/1E3)
    n_rejected <- 0
    list_index <- which(indexOut != 0, arr.ind = TRUE)
    for(i in seq_len(nrow(list_index))){
        row <- list_index[i,'row']
        col <- list_index[i,'col']
        fc <- zScore * datasd[row,col]
        clcount <- indexOut[row,col] * fc + lmu[row,col]
        
        #multiply size factor again
        art_out <- round(sf[col]*2^clcount)
        if(art_out < max_out){
            k[row,col] <- art_out
        }else{
            #remove super large outliers
            indexOut[row,col] <- 0 
            n_rejected <- n_rejected + 1
        }
    }
    mode(k) <- "integer"
    
    assay(ods, 'trueCounts', withDimnames=FALSE) <- counts(ods)
    counts(ods, withDimnames=FALSE) <- k
    assay(ods, "trueOutliers", withDimnames=FALSE) <- indexOut
    
    return(ods)
}

#' 
#' Aproximates standard deviation of counts in log2 space.
#' 
#' @noRd
sdLogScale <- function(mu, theta){
    sqrt(mu*(1 + mu/theta)) / ((mu + 1)*log(2)) * 2
}
