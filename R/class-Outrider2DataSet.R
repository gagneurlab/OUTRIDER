#' 
#' Outrider2DataSet class and constructors 
#' 
#' The Outrider2DataSet class is designed to store the whole 
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
#' @return An Outrider2DataSet object.
#' 
#' @author Christian Mertes \email{mertes@@in.tum.de}, 
#'             Felix Brechtmann \email{brechtma@@in.tum.de}
#' 
#' @rdname Outrider2DataSet-class
#' 
#' @examples
#' 
#' ods <- makeExampleOutrider2DataSet()
#' ods
#' 
#' ods <- makeExampleOutriderDataSet()
#' ods
#'    
setClass("Outrider2DataSet", contains="RangedSummarizedExperiment",
         slots = list(
             distribution    = "character",
             preprocessing   = "character",
             transformation  = "character",
             sf_norm = "logical"
         ),
         prototype = list(
             distribution    = "gaussian",
             preprocessing   = "log",
             transformation  = "none",
             sf_norm = TRUE
         )
)

#' check sample annotation within the colData slot of the SE object
#' @noRd
validateAssays <- function(object) {
    if(length(assayNames(object)) > 1){
        if(any(duplicated(assayNames(object)))){
            return(
                "OUTRIDER enforces unique assay names! Please provide such names."
            )
        }
    }
    NULL
}

checkNames <- function(object){
    n <- rownames(object)
    if(is.null(n) | any(duplicated(n))){
        return("Please provide unique rownames.")
    }
    n <- colnames(object)
    if(is.null(n) | any(duplicated(n))){
        return("Please provide unique colnames.")
    }
    cold <- colData(object) 
    if(!"sampleID" %in% colnames(cold) | any(duplicated(cold[,"sampleID"]))){
        return("Please provide unique sampleIDs in colData.")
    }
}

validateModelParameters <- function(object){
    # check distribution parameter
    if(!isScalarCharacter(object@distribution)){
        return("The distribution parameter should be a scalar character.")
    }
    if(!(object@distribution %in% c("negative binomial", "gaussian"))){
        return("Distribution should be one of 'negative binomial' and  
                'gaussian'. Other distributions are not yet implemented.")
    }
    # check transformation parameter
    if(!isScalarCharacter(object@transformation)){
        return("The transformation parameter should be a scalar character.")
    }
    if(!(object@transformation %in% c("log", "none"))){
        return("Transformation should be one of 'log' or 'none'. 
               Other transformations are not yet implemented.")
    }
    # check preprocessing parameter
    if(!isScalarCharacter(object@preprocessing)){
        return("The preprocessing parameter should be a scalar character.")
    }
    if(!(object@preprocessing %in% c("log", "vst", "none"))){
        return("Preprocessing should be one of 'log', 'vst' and 'none'. 
               Other preprocessing schemes are not yet implemented.")
    }
    # check sizefactor normalization yes/no
    if(!isScalarLogical(object@sf_norm)){
        return("The sf_norm parameter should be a scalar logical.")
    }
    NULL
}


#' general validate function
#' @noRd
validateOutrider2DataSet <- function(object) {
    c(
        checkNames(object),
        validateAssays(object),
        validateModelParameters(object)
    )
}
setValidity("Outrider2DataSet", validateOutrider2DataSet)


#' show method for ModelDataSet
#' @noRd
showOutrider2DataSet <- function(object) {
    cat("class: Outrider2DataSet\n")
    show(as(object, "RangedSummarizedExperiment"))
    cat("------------------- Model parameters -------------------\n")
    cat(paste0("Distribution:               ", 
               modelParams(object, "distribution")), "\n")
    cat(paste0("Sizefactor normalization:   ", 
               modelParams(object, "sf_norm")), "\n")
    cat(paste0("Preprocessing:              ", 
               modelParams(object, "preprocessing")), "\n")
    cat(paste0("Transformation:             ", 
               modelParams(object, "transformation")), "\n")
    cat("\n")
}

setMethod("show", "Outrider2DataSet", function(object) {
    showOutrider2DataSet(object)
})

#' @rdname Outrider2DataSet-class
#' @export
Outrider2DataSet <- function(se, inputData, colData, 
            distribution=c("gaussian", "negative binomial"), 
            sf_norm=TRUE, preprocessing=c("log", "vst", "none"), 
            transformation=c("none", "log"), ...){ 
    
    # check input
    distribution   <- tolower(distribution)
    distribution   <- match.arg(distribution)
    transformation <- tolower(transformation)
    transformation <- match.arg(transformation)
    preprocessing  <- tolower(preprocessing)
    preprocessing  <- match.arg(preprocessing)
    if(!isScalarLogical(sf_norm)){
        stop("sf_norm should be either TRUE or FALSE")
    }
    
    # use SummarizedExperiment object
    if(!missing(se)){
        if(!is(se, "SummarizedExperiment")){
            stop("'se' must be a SummarizedExperiment object")
        }
        
        # use raw count data
    } else if(!missing(inputData)){
        if(missing(colData)){
            cols <- colnames(inputData)
            if(is.null(cols)){
                cols <- paste0("sample", seq_len(ncol(inputData)))
            }
            colData <- DataFrame(sampleID=cols)
            rownames(colData) <- colData[['sampleID']]
            colnames(inputData) <- colData[['sampleID']]
        }
        inputData <- as.matrix(inputData)
        se <- SummarizedExperiment(
            assays=SimpleList(observed=inputData),
            colData=colData, 
            ...)
        
        # nothing provided
    } else {
        stop("At least one of the se or inputData argument has to be defined!")
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
    
    obj <- new("Outrider2DataSet", se, 
               distribution=distribution, transformation=transformation, 
               preprocessing=preprocessing, sf_norm=sf_norm)
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
    ods <- Outrider2DataSet(inputData=intensityData, distribution="gaussian", 
                            sf_norm=TRUE, preprocessing="log", 
                            transformation="none")
    
    assay(ods, "trueMeanLog", withDimnames=FALSE) <- mu
    assay(ods, "trueSdLog", withDimnames=FALSE) <- matrix(logSd, nrow=n, ncol=m)
    colData(ods)[['trueSizeFactor']]         <- sf
    metadata(ods)[['optimalEncDim']]         <- q
    metadata(ods)[['encDimTable']]           <- data.table(
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
    datasd <- assay(ods, "trueSdLog")
    lmu <- assay(ods, "trueMeanLog")
    
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
    
    assay(ods, 'trueObservations', withDimnames=FALSE) <- assay(ods, 
                                                               "observed")
    assay(ods, 'observed', withDimnames=FALSE) <- exp(k)
    assay(ods, "trueOutliers", withDimnames=FALSE) <- indexOut
    
    return(ods)
}


