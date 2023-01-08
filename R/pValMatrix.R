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
#' @param subsets A named list of named lists specifying any number of gene 
#'             subsets (can differ per sample). For each subset, FDR correction
#'             will be limited to genes in the subset, and the FDR corrected 
#'             pvalues stored as an assay in the ods object in addition to the 
#'             transcriptome-wide FDR corrected pvalues. See the examples for 
#'             how to use this argument. 
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
#' # example of restricting FDR correction to subsets of genes of interest
#' genesOfInterest <- list("sample_1"=sample(rownames(ods), 3), 
#'                          "sample_2"=sample(rownames(ods), 8), 
#'                          "sample_6"=sample(rownames(ods), 5))
#' ods <- computePvalues(ods, subsets=list("exampleSubset"=genesOfInterest))
#' padj(ods, subsetName="exampleSubset")[1:10,1:10]
#' ods <- computePvalues(ods, 
#'                  subsets=list("anotherExampleSubset"=rownames(ods)[1:5]))
#' padj(ods, subsetName="anotherExampleSubset")[1:10,1:10]
#' 
#' @exportMethod computePvalues
setGeneric("computePvalues", 
        function(object, ...) standardGeneric("computePvalues"))

#' @rdname computePvalues
#' @export
setMethod("computePvalues", "OutriderDataSet", function(object, 
            alternative=c("two.sided", "greater", "less"), method='BY', 
            subsets=NULL, BPPARAM=bpparam()){
    
    alternative <- match.arg(alternative)
    object <- pValMatrix(object, alternative, BPPARAM=BPPARAM)
    
    if(method != 'None'){
        object <- padjMatrix(object, method)
        
        # calculate FDR for each provided subset and assign to ods
        if(!is.null(subsets)){
            stopifnot(is.list(subsets))
            stopifnot(!is.null(names(subsets)))
            for(setName in names(subsets)){
                geneListSubset <- subsets[[setName]]
                object <- padjOnSubset(ods=object, method=method, 
                                    BPPARAM=BPPARAM,
                                    genesToTest=geneListSubset,
                                    subsetName=setName)
            }
        }
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
#'     will be stored as an assay with name 'padjust_'+subsetName.
#' @return OutriderDataSet containing adjusted pValues
#' 
#' @noRd
padjOnSubset <- function(ods, genesToTest, subsetName, method='BY', 
                            BPPARAM=bpparam()){
    
    # check input
    stopifnot(!is.null(genesToTest))
    stopifnot(is.list(genesToTest) || is.vector(genesToTest))
    
    # replicate subset genes for each sample if given as vector
    if(!is.list(genesToTest)){
        genesToTest <- rep(list(genesToTest), ncol(ods))
        names(genesToTest) <- colnames(ods)
    }
    
    # check that names are present and correspond to samples in ods
    stopifnot(!is.null(names(genesToTest)))
    if(!all(names(genesToTest) %in% colnames(ods))){
        stop("names(genesToTest) need to be sampleIDs in the given ods object.")
    }
    
    # compute FDR on the given subsets of genes
    fdrSubset <- bpmapply(colnames(ods), FUN=function(sampleId){
                                         
            # get genes to test for this sample
            genesToTestSample <- genesToTest[[sampleId]]
            padj <- rep(NA, nrow(ods))
            
            # if no genes present in the subset for this sample, return NAs
            if(is.null(genesToTestSample)){
                return(padj)
            }
            
            # get idx of junctions corresponding to genes to test
            if(is.character(genesToTestSample)){
                rowIdx <- sort(which(rownames(ods) %in% genesToTestSample))
            } else{
                stop("Genes in the list to test must be a character vector ",
                    "of geneIDs.")
            }
                             
            # check that rowIdx is not empty vector
            if(length(rowIdx) == 0){
                warning("No genes from the given subset found in the ods ",
                            "object for sample: ", sampleId)
                return(padj)
            }
                                         
            # retrieve pvalues of genes to test
            p <- pValue(ods)[rowIdx, sampleId]
                     
            # FDR correction on subset
            padjSub <- p.adjust(p, method=method)
                                         
            # return FDR on subset, filled with NA for all other genes
            padj[rowIdx] <- padjSub
            return(padj)
            
        }, SIMPLIFY=TRUE, BPPARAM=BPPARAM)
    rownames(fdrSubset) <- rownames(ods)
    colnames(fdrSubset) <- colnames(ods)
    
    # add FDR subset info to ods object and return
    padj(ods, subsetName=subsetName) <- fdrSubset
    validObject(ods)
    return(ods)
}

