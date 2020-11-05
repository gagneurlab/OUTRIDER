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
             fitModel        = "character"
         ),
         prototype = list(
             distribution    = "Gaussian",
             preprocessing   = "None",
             transformation  = "sf-log",
             fitModel        = "Autoencoder"
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
    if(!(object@distribution %in% c("Negative-Binomial", "Gaussian"))){
        return("Distribution should be one of Negative-Binomial and Gaussian. 
               Other distributions are not yet implemented.")
    }
    # check transformation parameter
    if(!isScalarCharacter(object@transformation)){
        return("The transformation parameter should be a scalar character.")
    }
    if(!(object@transformation %in% c("sf-log", "log", "None"))){
        return("Transformation should be one of sf-log, log or None. 
               Other transformations are not yet implemented.")
    }
    # check preprocessing parameter
    if(!isScalarCharacter(object@preprocessing)){
        return("The preprocessing parameter should be a scalar character.")
    }
    if(!(object@preprocessing %in% c("sf-log", "None"))){
        return("Preprocessing should be one of sf-log and None. 
               Other preprocessing schemes are not yet implemented.")
    }
    # check fit parameter
    if(!isScalarCharacter(object@fitModel)){
        return("The fitModel parameter should be a scalar character.")
    }
    if(!(object@fitModel %in% c("Autoencoder", "PCA"))){
        return("fitModel should be one of Negative-Binomial and Gaussian. 
               Other fitModels are not yet implemented.")
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
    cat(paste0("Preprocessing:              ", 
               modelParams(object, "preprocessing")), "\n")
    cat(paste0("Transformation:             ", 
               modelParams(object, "transformation")), "\n")
    cat(paste0("Fit-model:                  ", 
               modelParams(object, "fitModel"), "'"), "\n")
    cat("\n")
}

setMethod("show", "Outrider2DataSet", function(object) {
    showOutrider2DataSet(object)
})

#' @rdname OutriderModelDataSet-class
#' @export
Outrider2DataSet <- function(se, inputData, colData, ...) { # distr=NULL, trans=NULL, prepro="None", fitModel="Autoencoder",
    
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
            assays=SimpleList(raw=inputData),
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
    
    obj <- new("Outrider2DataSet", se) # distribution=distr, transformation=trans, preprocessing=prepro, fitModel=fitModel
    metadata(obj)[["version"]] <- packageVersion("OUTRIDER")
    validObject(obj)
    
    return(obj)
}


