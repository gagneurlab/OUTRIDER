#' 
#' Calculate P-values
#' 
#' This function computes the P-values based on fitted negative binomial model.
#' It computes two matrices with the same dimension as the count matrix
#' (samples x genes), which contain the corresponding P-values and 
#' adjusted P-values to every count. 
#' 
#' @param object An OutriderDataSet
#' @param alternative Can be one of "two.sided", "greater" or "less" to specify
#'             the alternative hypothesis used to calculate the P-values,
#'             defaults to "two.sided"
#' @param method Method used for multiple testing
#' @param modelFile The file name where the fit model should be taken from
#' @param BPPARAM Can be used to parallelize the computation, defaults to
#'             bpparam()
#' @param ... additional params, currently not used.
#' 
#' @return An OutriderDataSet object with computed nominal and adjusted P-values
#' 
#' @seealso p.adjust
#' 
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
                    method='BY', modelFile=NULL, BPPARAM=bpparam()){
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
    ctsData <- counts(ods)
    disp    <- mcols(ods)[,"disp"]
    mu      <- mcols(ods)[,"mu"]
    normF   <- normalizationFactors(ods)
    if(is.null(normF)){
        normF <- sizeFactors(ods)
    }
    
    pValMat <- bplapply(seq_along(ods), pValNorm, ctsData=ctsData, disp=disp,
            mu=mu, normF=normF, alternative=alternative, BPPARAM=BPPARAM)
    
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
    pless <- pnbinom(x, size = size, mu = s*mu)
    if(alternative == "less"){
        return(pless)
    }
    
    dval <- dnbinom(x, size = size, mu = s*mu)
    if(alternative ==  "greater"){
        return(1 - pless + dval)
    }
    
    return(2 * vapply(seq_along(pless), function(i) {
            min(0.5, pless[i], 1 - pless[i] + dval[i]) }, numeric(1)))
}


#' pValues per gene
#' computes pValues of all samples for a given gene.
#' 
#' @param index of gene
#' @param ods object
#'
#' @return vector of pValues
#' @noRd
pValNorm <- function(index, ctsData, disp, mu, normFact, alternative){
    data <- ctsData[index,]
    disp <- disp[index]
    mu   <- mu[index]
    if(is.matrix(normFact)){
        normFact <- normFact[index,]
    }
    pVal(data, size=disp, mu=mu, normFact, alternative=alternative)
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
