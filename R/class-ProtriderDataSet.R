#' 
#' ProtriderDataSet class and constructors 
#' 
#' The ProtriderDataSet class is designed to store the whole 
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
#' @return An ProtriderDataSet object.
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}, 
#'             Felix Brechtmann \email{brechtma@@in.tum.de}
#' 
#' @rdname ProtriderDataSet-class
#' 
#' @examples
#' 
#' ods <- makeExampleProtriderDataSet()
#' ods
#' 
#' ods <- makeExampleProtriderDataSet(dataset="Kremer")
#' ods
#'    
setClass("ProtriderDataSet", contains="Outrider2DataSet", 
         prototype = list(
             distribution    = "Gaussian",
             preprocessing   = "sf-log",
             transformation  = "None",
             fitModel        = "Autoencoder"
         ))

#' general validate function
#' @noRd
validateProtriderDataSet <- function(object) {
    c(
        checkNames(object),
        validateAssays(object),
        validateModelParameters(object)
    )
}
setValidity("ProtriderDataSet", validateProtriderDataSet)


#' show method for ProtriderDataSet
#' @noRd
showProtriderDataSet <- function(object) {
    cat("class: ProtriderDataSet\n")
    show(as(object, "Outrider2DataSet"))
}

setMethod("show", "ProtriderDataSet", function(object) {
    showProtriderDataSet(object)
})

#' @rdname ProtriderDataSet-class
#' @export
ProtriderDataSet <- function(se, intensityData, colData, ...) {
    
    # use SummarizedExperiment object
    if(!missing(se)){
        if(!is(se, "SummarizedExperiment")){
            stop("'se' must be a SummarizedExperiment object")
        }
        
        # use raw count data
    } else if(!missing(intensityData)){
        if(missing(colData)){
            cols <- colnames(intensityData)
            if(is.null(cols)){
                cols <- paste0("sample", seq_len(ncol(intensityData)))
            }
            colData <- DataFrame(sampleID=cols)
            rownames(colData) <- colData[['sampleID']]
            colnames(intensityData) <- colData[['sampleID']]
        }
        intensityData <- as.matrix(intensityData)
        se <- SummarizedExperiment(
            assays=SimpleList(intensities=intensityData), 
            colData=colData, 
            ...)
        
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
    se <- as(se, "RangedSummarizedExperiment")
    
    obj <- new("ProtriderDataSet", se)
    metadata(obj)[["version"]] <- packageVersion("OUTRIDER")
    validObject(obj)
    
    return(obj)
}

#' 
#' Create example data sets for OUTRIDER with protein intensities
#' 
#' Creates an example data set by simulating a data set based 
#' on random intensities following a log-normal distribution with injected
#' outliers with a fixed z score away from the mean of the gene.
#' 
#' @param n Number of simulated proteins 
#' @param m Number of simulated samples
#' @param freq Frequency of in-silico outliers
#' @param zScore Absolute z score of in-silico outliers (default 6).
#' @param inj Determines whether outliers are injected with the strategy 
#'            ('both', 'low', 'high'), default is 'both'.
#' @param sf Artificial Size Factors
#' @param q number of simulated latend variables.
#' 
#' @return An ProtriderDataSet containing an example dataset
#' 
#' @examples
#' # A generic dataset 
#' pds <- makeExampleProtriderDataSet()
#' pds
#' 
#' @export
makeExampleProtriderDataSet <- function(n=200, m=80, q=10, freq=1E-3, zScore=6,
                                       inj=c('both', 'low', 'high'), 
                                       sf=rnorm(m, mean=1, sd=0.1)){
    # make example data set 
    logSd <- rlnorm(n, meanlog=log(0.1), sdlog=0.5)  # log sd
    logMean <- 1                                   # log offset for mean expres
    sdVec <- rep(0.5, m)                           # sd for H matrix
    
    #
    # Simulate covariates.
    #
    H_true <- matrix(rnorm(m*q), nrow=m, ncol=q)
    D_true <- matrix(rnorm(n*q, sd=sdVec), nrow=n, ncol=q)
    y_true <- D_true %*% t(cbind(H_true))
    mu     <- t(t(rlnorm(n, logMean, logSd) + y_true)*sf)
    
    #
    # Simulate intensity Matrix with specified means.
    #
    k <- matrix(rnorm(m*n, mean=mu, sd=logSd), nrow=n, ncol=m)
    
    #
    # Create Outrider data set
    #
    intensityData <- exp(k)
    row.names(intensityData) <- paste0("feature_", seq_len(n))
    colnames(intensityData) <- paste0("sample_", seq_len(m))
    pds <- ProtriderDataSet(intensityData=intensityData)
    
    assay(pds, "trueMeanLog", withDimnames=FALSE) <- mu
    assay(pds, "trueSdLog", withDimnames=FALSE) <- matrix(logSd, nrow=n, ncol=m)
    colData(pds)[['trueSizeFactor']]         <- sf
    metadata(pds)[['optimalEncDim']]         <- q
    metadata(pds)[['encDimTable']]           <- data.table(
        encodingDimension=q, evaluationLoss=1, evalMethod='simulation')
    
    #
    # inject outliers
    #
    indexOut <- matrix(nrow=n, ncol=m,
                        sample(c(-1,1,0), m*n, replace=TRUE, 
                        prob=c(freq/2, freq/2, 1-freq)))
    indexOut <- switch(match.arg(inj),
                        low  = -abs(indexOut),
                        high =  abs(indexOut),
                        indexOut
    )
    # normtable <- t(t(k)/sf)
    datasd <- assay(pds, "trueSdLog")
    lmu <- assay(pds, "trueMeanLog")
    
    # inject outliers
    max_out <- 1E2 * min(max(k), .Machine$double.xmax/1E3)
    n_rejected <- 0
    list_index <- which(indexOut != 0, arr.ind = TRUE)
    for(i in seq_len(nrow(list_index))){
        row <- list_index[i,'row']
        col <- list_index[i,'col']
        fc <- zScore * datasd[row,col]
        corK <- indexOut[row,col] * fc + lmu[row,col]
        
        #multiply size factor again
        art_out <- sf[col]*corK
        if(art_out < max_out){
            k[row,col] <- art_out
        }else{
            #remove super large outliers
            indexOut[row,col] <- 0 
            n_rejected <- n_rejected + 1
        }
    }
    
    assay(pds, 'trueIntensities', withDimnames=FALSE) <- assay(pds, 
                                                                "intensities")
    assay(pds, 'intensities', withDimnames=FALSE) <- exp(k)
    assay(pds, "trueOutliers", withDimnames=FALSE) <- indexOut
    
    return(pds)
}

