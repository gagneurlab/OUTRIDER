#'
#' Checks the count requirements so it has enough counts per gene and sample
#' 
#' @noRd
checkCountRequirements <- function(ods, test=FALSE){
    checkOutriderDataSet(ods)
    if(ncol(ods) == 0){
        stop("Please provide at least one sample.")
    }
    if(nrow(ods) == 0){
        stop("Please provide at least one gene.")
    }
    filterGenes <- rowSums(counts(ods)) == 0
    if(any(filterGenes) & isFALSE(test)){
        stop("There are genes without any read. Please filter first ",
                "the data with: ods <- filterExpression(ods)")
    }
    filterGenes <- filterGenes | rowSums(counts(ods)) < ncol(ods)/100
    if(any(filterGenes) & isFALSE(test)){
        stop("The model requires for each gene at least 1 read ", 
                "in every 100 sample. Please filter first the data with: ",
                "ods <- filterExpression(ods)")
    }
    return(invisible(filterGenes))
}

#' 
#' Checks that is is an OutriderDataSet object
#' @noRd
checkOutriderDataSet <- function(ods){
    if(!is(ods, 'OutriderDataSet')){
        stop('Please provide an OutriderDataSet.')
    }
    return(invisible(TRUE))
}

#' 
#' Checks that the object contains sizeFactors
#' @noRd
checkSizeFactors <- function(ods, funName=sys.calls()[[1]]){
    checkOutrider2DataSet(ods)
    if(is.null(sizeFactors(ods))){
        stop("Please calculate the size factors before calling the '", funName,
                "' function. Please do: ods <- estimateSizeFactors(ods)")
    }
    return(invisible(TRUE))
}

#' 
#' Checks if the final padjust and zscore values are present
#' @noRd
checkFullAnalysis <- function(ods, funName=sys.calls()[[1]]){
    checkSizeFactors(ods)
    if(!'padjust' %in% assayNames(ods)){
        stop("Please calculate the P-values before calling the '", funName,
                "' function. Please do: ods <- computePvalues(ods)")
    }
    if(!'zScore' %in% assayNames(ods)){
        stop("Please calculate the Z-scores before calling the '", funName,
                "' function. Please do: ods <- computeZscores(ods)")
    }
    return(invisible(TRUE))
}

#' 
#' Checks if the provided thetaRange is correct
#' @noRd
checkThetaRange <- function(thetaRange){
    if(length(thetaRange) != 2){
        stop('Please provide a range for thetaRange (eg: `c(0.1, 250)`).')
    }
    if(!all(is.finite(thetaRange))){
        stop('Please provide finite values for the theta range.')
    }
    if(thetaRange[1] > thetaRange[2]){
        stop('The first element of the range has to be smaller', 
                ' than the second one.')
    }
    return(invisible(TRUE))
}

#' 
#' Checks general data requirements (for general Outrider2DataSet)
#' @noRd
checkDataRequirements <- function(ods, test=FALSE){
    checkOutrider2DataSet(ods)
    if(ncol(ods) == 0){
        stop("Please provide at least one sample.")
    }
    if(nrow(ods) == 0){
        stop("Please provide at least one feature.")
    }
    filterFeatures <- rowSums(is.na(raw(ods))) > 0
    if(any(filterFeatures) & isFALSE(test)){
        stop("There are features that contain NAs. Please filter first ",
             "the data with: ods <- filterData(ods)")
    }
    return(invisible(filterFeatures))
}

checkPreprocessing <- function(ods){
    if(modelParams(ods, "preprocessing") != "None"){
        if(!("preprocessed" %in% assayNames(ods))){
            stop("The ods needs to be preprocessed first. Please run first: ", 
                 "ods <- preprocess(ods)")
        }
    }
    if(is(ods, "OutriderDataSet")){
        checkSizeFactors(ods, funName=sys.calls()[[1]])
    }
    return(invisible(TRUE))
}

#' 
#' Checks that is is an Outride2rDataSet object
#' @noRd
checkOutrider2DataSet <- function(ods){
    if(!is(ods, 'Outrider2DataSet')){
        stop('Please provide an Outrider2DataSet.')
    }
    return(invisible(TRUE))
}

#' 
#' Checks that is is an OutriderDataSet object
#' @noRd
checkFitInR <- function(ods, fitInR){
    checkOutrider2DataSet(ods)
    if(!is(ods, 'OutriderDataSet') & isTRUE(fitInR)){
        stop('The fitting function in R works only for OutriderDataSet. ",
             "For the more general Outrider2Dataset, please use the python "
             "version by setting usePython=TRUE.')
    }
    return(invisible(TRUE))
}
