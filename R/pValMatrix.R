#' 
#' This is the main function to calculate all P-values.
#' 
#' Computes two matrices with the same dimension as the count matrix, 
#' which contain the corresponding P-values and adjusted P-values to every 
#' count.
#' 
#' 
#' @param object object
#' @param alternative "two.sided", "greater" or "less" P-values can 
#'             be computed, default is "two.sided"
#' @param method method used for multiple testing
#' @param modelFile The file name where the fit model should be taken from
#' @param BPPARAM by default bpparam()
#' @param ... additional params, currently not used.
#' 
#' @return An OutriderDataSet object with computed P-values and padjusted values
#' 
#' @seealso p.adjust
#' @docType methods
#' @name computePvalues
#' @rdname computePvalues
#' @aliases computePvalues computePvalues,OutriderDataSet-method
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods <- estimateSizeFactors(ods)
#' ods <- fit(ods)
#' ods <- computePvalues(ods)
#' 
#' assays(ods)[['pValue']][1:10,1:10]
#' 
#' ods <- computeZscores(ods)
#' res <- results(ods)
#' res
#' 
#' @exportMethod computePvalues
setGeneric("computePvalues", 
        function(object, ...) standardGeneric("computePvalues"))

#' @rdname computePvalues
#' @export
setMethod("computePvalues", "OutriderDataSet", function(object, 
                    alternative=c("two.sided", "greater", "less"), 
                    method='BH', modelFile=NULL, BPPARAM=bpparam()){
    if(!is.null(modelFile)){
        object <- readNBModel(object, modelFile)
        object <- estimateSizeFactors(object)
    }
    alternative <- match.arg(alternative)
    object <- pValMatrix(object, alternative, BPPARAM=BPPARAM)
    padjMatrix(object, method)
})

pValMatrix <- function(ods, alternative, BPPARAM){ 
    if(!all(c("disp", "mu") %in% colnames(mcols(ods)))){
        stop(paste("Please fit the models first to estimate disp and",
                "mu by running:\n\tods <- fit(ods)"))
    }
    pValMat <- bplapply(1:length(ods), pValNorm, 
            ods=ods, alternative=alternative, BPPARAM=BPPARAM)
    pValMat <- matrix(unlist(pValMat), nrow=length(ods), byrow=TRUE)
    assays(ods)[["pValue"]] <- pValMat
    validObject(ods)
    return(ods)
}


#' 
#' Internal P-value calculation
#' 
#' Since the distribution has an integer support one needs to take care of 
#' summing over the x bin on both sides. Further need to take care, that 
#' P-values do not become larger then one.
#'
#' \deqn{
#' p_{ij} = 2 \cdot min \left\{
#'         \frac{1}{2},  
#'         \sum_{0}^{k_{ij}} NB(\mu_i\cdot c_{ij} ,\theta_i),
#'         1 - \sum_{0}^{k_{ij-1}} NB(\mu_i\cdot c_{ij} ,\theta_i)
#'     \right\}
#' }
#' 
#' @param x count 
#' @param size fitted dispersion (size) for the neg. binomial dist
#' @param mu fitted mean for the neg. binomial dist
#' @param s sizeFactor or normalizationFactor
#'
#' @return p Value
#' @noRd
pVal <- function(x, size, mu, s, alternative){
    if(alternative == "two.sided"){
        p <- 2 * min(0.5, pnbinom(x, size = size, mu = s*mu), 
            1 - pnbinom(x, size = size, mu = s*mu) + 
                dnbinom(x, size = size, mu = s*mu))
    }
    if(alternative == "less"){
        p <- pnbinom(x, size = size, mu = s*mu)
    }
    if(alternative ==  "greater"){
        p <- 1 - pnbinom(x, size = size, mu = s*mu) + 
            dnbinom(x, size = size, mu = s*mu)
    }
    return(p)
}


#' pValues per gene
#' computes pValues of all samples for a given gene.
#' 
#' @param index of gene
#' @param ods object
#'
#' @return vector of pValues
#' @noRd
pValNorm <- function(index, ods, alternative){
    data = counts(ods[index,])
    disp <- mcols(ods)[index,"disp"]
    mu <- mcols(ods)[index,"mu"]
    normF <- sizeFactors(ods)
    if(!is.null(normalizationFactors(ods))){
        normF <- normalizationFactors(ods[index,])
    }
    sapply(1:length(data), function(i) 
        pVal(data[i], size=disp, mu=mu, normF[i], alternative=alternative))
}


#' 
#' FDR correction
#' 
#' Function to correct the computed pValues using BY FDR correction.
#' This function is called by the computePvalues method.
#'
#' @param ods OutriderDataSet
#' @param method adjustment method, by default 'BH'
#' @return OutriderDataSet containing adjusted pValues
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' ods <- OUTRIDER(ods)
#' ods <- padjMatrix(ods, method='BH')
#' head(assays(ods)[["padjust"]][,1:10])
#'
#' ods <- padjMatrix(ods, method='fdr')
#' head(assays(ods[,1:10])[["padjust"]])
#' 
#' @export
padjMatrix <- function(ods, method='BH'){
    padj <- apply(assays(ods)[["pValue"]], 2, p.adjust, method=method)
    assays(ods)[["padjust"]] <- padj 
    validObject(ods)
    return(ods)
}
