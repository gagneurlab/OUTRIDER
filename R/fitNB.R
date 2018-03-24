#'
#' Fit a NB distribution to the genecounts per gene over all samples using the
#' precomputed control factors.
#' Parallized version.
#' 
#' @param object OutriderDataSet
#' @param modelFile file name to save the estimated parameters
#' @param BPPARAM by default bpparam()
#' @param ... additional arguments.
#' 
#' @docType methods
#' @name fit
#' @rdname fit
#' @aliases fit fit,OutriderDataSet-method
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
#' ods <- estimateSizeFactors(ods)
#' ods <- fit(ods)
#'   
#' mcols(ods)[1:10,c('mu', 'disp')]
#'   
#' @exportMethod fit 
setGeneric("fit", function(object, ...) standardGeneric("fit"))

#' @rdname fit
#' @export
setMethod("fit", "OutriderDataSet", function(object, dropExtremeRank=FALSE, 
                                             modelFile=NULL, BPPARAM=bpparam()){
    fitNB(object, dropExtremeRank=dropExtremeRank, modelFile=modelFile, 
          BPPARAM=BPPARAM)
})


fitNB <- function(ods, dropExtremeRank, modelFile, BPPARAM){
    fitparameters <- bplapply(1:length(ods), fitNegBinom, 
            ods=ods, dropExtremeRank=dropExtremeRank, BPPARAM=BPPARAM)
    fitparameters <- DataFrame(t(matrix(unlist(fitparameters), nrow = 2)))
    mcols(ods)[c('mu', 'disp')] <- fitparameters
    validObject(ods)
    
    # write out model file if requested
    if(!is.null(modelFile)){
        writeNBModel(ods, modelFile)
    }
    return(ods)
}


## functions called inside the fit

#initial value estimation. 
initialSizeMu <- function(data, norm){
    m <- mean(data/norm)
    v <- var(as.vector(data/norm)) # inserted as.vector
    size <- if(v > m){
        m^2/(v - m)
    } else {
        1  
    }  
    c("size"=c(size), "mu"=c(m))
}

# log of the probability denity function
dens <- function(x, size, mu, s) dnbinom(x, size=size, mu=mu*s, log=TRUE)

# log likelihood
loglikelihood <- function(sizemu, x, SizeF){
    -sum(sapply(1:length(x), function(j){
        dens(x[j], size=sizemu[1], mu=sizemu[2], s=SizeF[j]) }))
}

# gradient of log likelihood
gradloglikelihood <- function(sizemu, x, SizeF){
    r <- sizemu[1]
    m <- sizemu[2]
    s <- SizeF
    c(- sum(log(r/(m*s+r))) - sum((m*s-x)/(r+m*s)) 
      - sum(digamma(x+r)) + length(x) * digamma(r),
      -r/m * sum((x-m*s)/(r+m*s)))
}


#fit for individual gene
fitNegBinom <- function(index, ods, dropExtremeRank){
    data <- counts(ods[index,])[1,]
    
    normF <- sizeFactors(ods)
    if(!is.null(normalizationFactors(ods))){
        normF <- normalizationFactors(ods[index,])
    }
    if(is.null(normF)){
        stop(paste("Please provide sizeFactors or normalizationFactors", 
                "for better estimates!\n  To overcome this just set the", 
                "sizeFactors like this:\n    sizeFactors(ods) <- 1"))
    }
    if(dropExtremeRank == TRUE){
        #option to remove extreme points from distribution prior to fitting.
        r <- rank(data)
        ndrop <- round(length(r)*0.01)
        data <- data[r>ndrop & r<length(r)-ndrop]
        normF <- normF[r>ndrop & r<length(r)-ndrop]
        # if(NF != 1){
        #     NF <- NF[r>ndrop & r<length(r)-ndrop]
        # }
    }
    
    ##correct s factor
    est <- optim(par=initialSizeMu(data, normF), fn=loglikelihood, 
        gr = gradloglikelihood, x=data, SizeF=normF, method="L-BFGS-B", 
        lower=c(0.01,0.01))
    c(est$par["mu"], est$par["size"])
} 





