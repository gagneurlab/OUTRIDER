#'
#' Fit the negative binomial distribution
#' 
#' Fit a negative binomial (NB) distribution to the counts per gene
#' over all samples using, if available, the precomputed control factors.
#' If no normalization factors are provided only the sizeFactors are used.
#' 
#' @param object An OutriderDataSet
#' @param BPPARAM by default bpparam()
#' @param ... additional arguments, currently not used.
#' @return An OutriderDataSet object with the fitted model. Accessible through:
#'             \code{mcols(ods)[,c('mu', 'theta')]}.
#' 
#' @docType methods
#' @name fit
#' @rdname fit
#' @aliases fit fit,OutriderDataSet-method
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- fit(ods)
#' 
#' mcols(ods)[1:10,c('mu', 'theta')]
#' 
#' @exportMethod fit 
setGeneric("fit", function(object, ...) standardGeneric("fit"))

#' @rdname fit
#' @export
setMethod("fit", "OutriderDataSet", function(object, BPPARAM=bpparam()){
    fitTheta(object, BPPARAM=BPPARAM)})

fitTheta <- function(ods, BPPARAM){
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    
    ctsData <- counts(ods)
    normF <- normalizationFactors(ods)
    if(is.null(normF)){
        normF <- sizeFactors(ods)
    }
    if(is.null(normF)){
        stop("Please provide sizeFactors or normalizationFactors for better ",
                "estimates!\n  To compute the sizeFactors please run at least:",
                " ods <- estimateSizeFactors(ods).")
    }
    
    excludeMask <- sampleExclusionMask(ods, aeMatrix=TRUE)
    
    fitparameters <- bplapply(seq_along(ods), fitNegBinom, normF=normF,
            ctsData=ctsData, excludeMask=excludeMask, BPPARAM=BPPARAM)
    
    mcols(ods)['mu'] <- vapply(fitparameters, "[[", double(1), "mu")
    theta(ods) <- vapply(fitparameters, "[[", double(1), "size")
    
    assay(ods)

    validObject(ods)
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

# log likelihood
loglikelihood <- function(sizemu, x, SizeF){
    -sum(dnbinom(x, size=sizemu[1], mu=sizemu[2]*SizeF, log=TRUE))
}

# gradient of log likelihood
gradloglikelihood <- function(sizemu, x, SizeF){
    r <- sizemu[1]
    m <- sizemu[2]
    s <- SizeF
    c(gradTheta(r, x, m*s),
        -r/m * sum((x-m*s)/(r+m*s)))
}

gradTheta <- function(theta, k, mu){
    ll <- log(mu+theta)-log(theta)-1 + (k+theta)/(theta+mu) - 
        digamma(k+theta) + digamma(theta)
    sum(ll)
}

#fit for individual gene
fitNegBinom <- function(index, ctsData, normF, excludeMask){
    data <- ctsData[index,]
    if(is.matrix(normF)){
        normF <- normF[index,]
    }
    stopifnot(!is.null(normF))
    
    data  <- data[ excludeMask[index,] == 1]
    normF <- normF[excludeMask[index,] == 1]
    
    ##correct s factor
    par <- initialSizeMu(data, normF)
    est <- list(par=c(mu=NA_real_, size=NA_real_))
    # try 3 times to find the optimum and change init values slightly each time
    for(i in seq_len(3)){
        est <- tryCatch(
                optim(par=par, fn=loglikelihood, gr=gradloglikelihood, x=data,
                        SizeF=normF, method="L-BFGS-B", lower=c(0.01,0.01)),
                error = function(e){
                        warning('Fit (', index, ') resulted in error using ', 
                                'init values: ', e$message)})
        if(!is.na(est$par["mu"])){
            break
        }
        par <- abs(par + rnorm(2, 0.2))
    }
    if(is.na(est$par["mu"])){
        # this is currently only happening on windows 32 bit. Probably due to
        # the lower precision compared to 64 bit.
        warning("Fit (", index, ") returned NA's. Using init values instead.")
        return(par)
    }
    c(est$par["mu"], est$par["size"])
}
