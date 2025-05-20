#' 
#' Autoencoder function to correct for confounders.
#' 
#' This is the wrapper function for the autoencoder implementation. 
#' It can be used to call the standard R implementation or the experimental
#' Python implementation.
#'
#' @param ods An OutriderDataSet object
#' @param q The encoding dimensions
#' @param implementation "autoencoder", the default, will use the autoencoder
#'             implementation. Also 'pca' and 'peer' can be used to control
#'             for confounding effects
#' @param BPPARAM A 
#'     \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'             instance to be used for parallel computing.
#' @param ... Further arguments passed on to the specific implementation method.
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' implementation <- 'autoencoder'
#' \dontshow{
#'     ods <- ods[1:10,1:10]
#'     implementation <- 'pca'
#' }
#' ods <- estimateSizeFactors(ods)
#' ods <- controlForConfounders(ods, implementation=implementation)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
#' @export
controlForConfounders <- function(ods, q,
                    implementation=c('autoencoder', 'pca'),
                    BPPARAM=bpparam(), ...){
    
    # error checking
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    checkSizeFactors(ods)
    
    if(!missing(q)){
        if(!(is.numeric(q) & q > 1 & q <= min(dim(ods)))){
            stop("Please provide for q an integer greater than 1 and smaller ", 
                    "than number of samples or genes.")
        }
    } else {
        q <- getBestQ(ods)
        if(is.na(q)){
            ods <- estimateBestQ(ods)
            q <- getBestQ(ods)
            message('Using estimated q with: ', q)
        } else {
            message('Using provided q with: ', q)
        }
    }
    
    # pass on to the approriate implementation
    implementation <- tolower(implementation)
    implementation <- match.arg(implementation)
    aeFun <- switch(implementation,
        autoencoder = { 
                function(ods, q, ...){ fitAutoencoder(ods, q, ...) } },
        pca         = { 
                function(ods, q, BPPARAM, ...){ autoCorrectPCA(ods, q, ...) } },
        stop("Requested control implementation is unknown.")
    )
    
    message(date(), ": Using the ", implementation, 
            " implementation for controlling.")
    ans <- aeFun(ods=ods, q=q, BPPARAM=BPPARAM, ...)
    message(date(), ": Used the ", implementation, 
            " implementation for controlling.")
    
    return(ans)
}

#' 
#' Extracting the latent space
#' 
#' Extracts the latent space from the OutriderDataSet object 
#' determined by the autoencoder.
#'
#' @param ods An OutriderDataSet
#'
#' @return A matrix containing the latent space determined by the autoencoder.
#'
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' \dontshow{
#'     ods <- ods[1:10, 1:10]
#' }
#' ods <- estimateSizeFactors(ods)
#' ods <- controlForConfounders(ods, implementation="pca")
#' computeLatentSpace(ods)[,1:6]
#' 
#' @export
computeLatentSpace <- function(ods){
    stopifnot(is(ods, 'OutriderDataSet'))
    if(is.null(D(ods)) | is.null(E(ods))){
        stop('The D or E weights are not computed yet. Please fit the ',
                'autoencoder before extracting the latent space.')
    }
    if(any(metadata(ods)[['dim']] != dim(ods))){
        stop('The OutriderDataSet dimension changed and does not match with ',
                'the existing autoencoder fit. Please refit the autoencoder.')
    }
    
    return(H(ods))
}