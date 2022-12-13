#' 
#' Calculate P-values
#' 
#' This function computes the P-values based on the fitted negative binomial
#' model. It computes two matrices with the same dimension as the count matrix
#' (samples x genes), which contain the corresponding P-values and 
#' adjusted P-values of every count. 
#' 
#' @param object An OutriderDataSet
#' @param alternative Can be one of "two.sided", "greater" or "less" to specify
#'             the alternative hypothesis used to calculate the P-values,
#'             defaults to "two.sided"
#' @param method Method used for multiple testing
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

#' @rdname computePvalues
#' @export
setMethod("computePvalues", "OutriderDataSet", function(object, 
            alternative=c("two.sided", "greater", "less"), method='BY', 
            BPPARAM=bpparam()){
    
    alternative <- match.arg(alternative)
    object <- pValMatrix(object, alternative, BPPARAM=BPPARAM)
    
    if(method != 'None'){
        object <- padjMatrix(object, method)
    }
    object
})

pValMatrix <- function(ods, alternative, BPPARAM){ 
    if(!all(c("theta") %in% colnames(mcols(ods)))){
        stop(paste("Please fit the models first to estimate theta and",
                "mu by running:\n\tods <- fit(ods)"))
    }
    ctsData <- counts(ods)
    mu <- rep(1, nrow(ods))
    if('mu' %in% names(mcols(ods))){
        mu <- mcols(ods)[,"mu"]
    }
    normF <- normalizationFactors(ods)
    if(is.null(normF)){
        normF <- sizeFactors(ods)
    }
    thetaMat <- outer(theta(ods), thetaCorrection(ods))
    
    pValMat <- bplapply(seq_along(ods), pVal, ctsData=ctsData, theta=thetaMat,
            mu=mu, normF=normF, alternative=alternative, BPPARAM=BPPARAM)
    
    pValue(ods) <- matrix(unlist(pValMat), nrow=length(ods), 
            byrow=TRUE, dimnames=dimnames(ods))
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
#' @param size fitted theta (size) for the neg. binomial dist
#' @param mu fitted mean for the neg. binomial dist
#' @param s sizeFactor or normalizationFactor
#'
#' @return p Value
#' @noRd
pVal <- function(index, ctsData, theta, mu, normFact, alternative){
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

#' 
#' FDR correction on subset of genes
#' 
#' Function to correct the computed pValues using BY FDR correction, limited to 
#' a subset of genes.
#'
#' @param ods OutriderDataSet
#' @param method adjustment method, by default 'BY'
#' @param genesToTest A named list with the subset of genes to test per sample.
#'     The names must correspond to the sampleIDs in the ods object.
#' @param subsetName The name under which the resulting FDR corrected pvalues 
#'     will be stored in metadata(ods) (after prefix 'FDR_').
#' @return OutriderDataSet containing adjusted pValues
#' 
#' @noRd
padjOnSubset <- function(ods, genesToTest, subsetName, method='BY', 
                            BPPARAM=bpparam()){
    
    # check input
    stopifnot(!is.null(genesToTest))
    stopifnot(is.list(genesToTest))
    stopifnot(!is.null(names(genesToTest)))
    if(!all(names(genesToTest) %in% colnames(ods))){
        stop("names(genesToTest) need to be sampleIDs in the given ods object.")
    }
    
    # compute FDR on the given subsets of genes
    FDR_subset <- rbindlist(bpmapply(names(genesToTest), genesToTest, 
        FUN=function(sample_id, genes_to_test_sample){
                                         
            # get idx of junctions corresponding to genes to test
            if(is.character(genes_to_test_sample)){
                row_idx <- sort(which(rownames(ods) %in% genes_to_test_sample))
            } else{
                stop("Genes in the list to test must be a character vector ",
                    "of geneIDs.")
            }
                             
            # check that row_idx is not empty vector
            if(length(row_idx) == 0){
                warning("No genes from the given subset found in the ods ",
                            "object for sample: ", sample_id)
                return(data.table(gene=character(0), 
                                    sampleID=character(0), 
                                    pval=numeric(0),
                                    FDR_subset=numeric(0), 
                                    idx=integer(0)))
            }
                                         
            # retrieve pvalues of junctions to test
            p <- pValue(ods)[row_idx, sample_id]
                     
            # FDR correction
            pa <- p.adjust(p, method=method)
                                         
            # return FDR on subset
            dt <- data.table(gene=rownames(ods)[row_idx], sampleID=sample_id, 
                                pval=p, FDR_subset=pa, 
                                idx=row_idx)
            return(dt)
        }, SIMPLIFY=FALSE, BPPARAM=BPPARAM))
    
    # add FDR subset info to ods object and return
    metadata(ods)[[paste("FDR", subsetName, sep="_")]] <- FDR_subset
    validObject(ods)
    return(ods)
}

