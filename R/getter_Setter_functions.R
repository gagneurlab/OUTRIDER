zScore <- function(ods){
    if(!'zScore' %in% assayNames(ods)){
        stop('Please compute first the Z-scores before retrieving them.')
    }
    assay(ods, 'zScore')
}

pValue <- function(ods){
    if(!'pValue' %in% assayNames(ods)){
        stop('Please compute first the P-values before retrieving them.')
    }
    assay(ods, 'pValue')
}

padj <- function(ods){
    if(!'padjust' %in% assayNames(ods)){
        stop('Please compute first the P-values before retrieving', '
                the adjusted ones.')
    }
    assay(ods, 'padjust')
}

#' @export dispersions
setMethod("dispersions", signature(object="OutriderDataSet"), function(object, ...){
    mcols(object)[['disp']]
})
