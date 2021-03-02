#' 
#' Calculate P-values
#' 
#' This function computes the P-values based on the fitted negative binomial
#' model. It computes two matrices with the same dimension as the count matrix
#' (samples x genes), which contain the corresponding P-values and 
#' adjusted P-values of every count. 
#' 
#' @param object An OutriderDataSet
#' @param distribution The distribution with which p values will be computed.
#'             Either 'nb' (for negative binomial) or 'gaussian'.
#' @param alternative Can be one of "two.sided", "greater" or "less" to specify
#'             the alternative hypothesis used to calculate the P-values,
#'             defaults to "two.sided"
#' @param method Method used for multiple testing. Set to \code{NULL} to 
#'             not run multiple testing correction.
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
#' ods <- makeExampleOutriderDataSet()
#' \dontshow{
#'     ods <- ods[1:10,1:10]
#' }
#' ods <- estimateSizeFactors(ods)
#' ods <- fit(ods)
#' 
#' ods <- computePvalues(ods)
#' 
#' assays(ods)[['pValue']][1:10,1:10]
#' 
#' @exportMethod computePvalues
setGeneric("computePvalues", 
        function(object, ...) standardGeneric("computePvalues"))

pValMatrix <- function(object, 
                        distribution=c("nb", "gaussian"), 
                        alternative=c("two.sided", "greater", "less"), 
                        method='BY', BPPARAM=bpparam()){
    
    alternative <- match.arg(alternative)
    distribution <- tolower(distribution)
    distribution <- match.arg(distribution)
    
    if(distribution == "nb"){
        object <- pValMatrix.nb(object, alternative, BPPARAM=BPPARAM)
    } else if(distribution == "gaussian"){
        object <- pValMatrix.gaussian(object, alternative, BPPARAM=BPPARAM)
    } else{
        stop("P value calculation for distribution ", distribution,
            " not implemented.")
    }
    metadata(object)$pvalDistribution <- distribution
    
    if(!is.null(method)){
        object <- padjMatrix(object, method)
    }
    object
}

#' @rdname computePvalues
#' @export
setMethod("computePvalues", "Outrider2DataSet", pValMatrix)

#' @rdname computePvalues
#' @export
setMethod("computePvalues", "OutriderDataSet", pValMatrix)

pValMatrix.nb <- function(ods, alternative, BPPARAM){ 
    if(!all(c("theta") %in% colnames(mcols(ods)))){
        stop(paste("Please fit the models first to estimate theta and",
                "mu by running:\n\tods <- fit(ods)"))
    }
    ctsData <- preprocessed(ods) # preprocessed == observed counts for outrider
    mu <- rep(1, nrow(ods))
    if('mu' %in% names(mcols(ods))){
        mu <- mcols(ods)[,"mu"]
    }
    normF <- normalizationFactors(ods)
    if(is.null(normF)){
        normF <- sizeFactors(ods)
    }
    thetaMat <- outer(theta(ods), thetaCorrection(ods))
    
    pValMat <- bplapply(seq_along(ods), pVal.nb, ctsData=ctsData, 
            theta=thetaMat, mu=mu, normF=normF, alternative=alternative, 
            BPPARAM=BPPARAM)
    
    pValue(ods) <- matrix(unlist(pValMat), nrow=length(ods), 
            byrow=TRUE, dimnames=dimnames(ods))
    validObject(ods)
    return(ods)
}


#' 
#' Internal NB P-value calculation
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
#' @param size fitted theta (size) for the neg. binomial dist
#' @param mu fitted mean for the neg. binomial dist
#' @param s sizeFactor or normalizationFactor
#'
#' @return p Value
#' @noRd
pVal.nb <- function(index, ctsData, theta, mu, normFact, alternative){
    x    <- ctsData[index,]
    mu   <- mu[index]
    
    # check if you did not get only sizeFactors
    if(is.matrix(normFact)){
        normFact <- normFact[index,]
    }
    size <- theta[index,]
    
    pless <- pnbinom(x, size=size, mu=normFact*mu)
    if(alternative == "less"){
        return(pless)
    }
    
    dval <- dnbinom(x, size=size, mu=normFact*mu)
    if(alternative ==  "greater"){
        return(1 - pless + dval)
    }
    
    return(2 * pmin(0.5, pless, 1 - pless + dval))
}

#' Gaussian P value calculation
#' @noRd
pValMatrix.gaussian <- function(ods, alternative, BPPARAM){ 
    
    inputData <- preprocessed(ods)
    predicted <- normalizationFactors(ods)
    residuals <- inputData - predicted
    
    sd <- rowSds(residuals, na.rm=TRUE)
    pval <- pnorm(inputData, mean=predicted, 
        sd=matrix(sd, nrow=nrow(inputData), ncol=ncol(inputData), byrow=FALSE))
    
    if(alternative == "less"){
        pValue(ods) <- pval
    }
    if(alternative == "greater"){
        pValue(ods) <- 1 - pval
    } 
    else{
        pval <- 2 * pmin(pval, 1 - pval)
        pValue(ods) <- pval
        
    }
    
    validObject(ods)
    return(ods)
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
#' @noRd
padjMatrix <- function(ods, method='BH'){
    padj(ods) <- apply(pValue(ods), 2, p.adjust, method=method)
    
    validObject(ods)
    return(ods)
}

