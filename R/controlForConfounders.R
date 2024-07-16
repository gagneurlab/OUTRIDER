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
            q <- estimateBestQ(ods)
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

#'
#' Estimation of Q
#' 
#' Estimating the best q for the given data set by Optimal Hard Thresholding
#' 
#' @param ods An OutriderDataSet object
#' @param zScores A z-score matrix
#' @import RMTstat
#' @import pracma
#' @return The estimated dimension of hidden confounders
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' 
#' estimateBestQ(ods)
#' 
#' @export
estimateBestQ <- function(ods=NULL, zScores=NULL){
  if (is.null(ods) & is.null(zScores)){
    stop("Please provide an OutriderDataSet or a z-score matrix.")
  }
  else if (!is.null(ods)){
    OUTRIDER:::checkOutriderDataSet(ods)
    if (!is.null(zScores)){
      warning("Provided z-scores are ignored and recalculated from ods.")
    }
    
    # Control for sequencing depth
    if (is.null(sizeFactors(ods))){
      ods <- estimateSizeFactors(ods)
    }
    controlledCounts <- t(t(counts(ods, normalized=FALSE)) / sizeFactors(ods)) 
    
    # Log-transform controlled counts
    logControlCounts <- log2((controlledCounts +1) / (rowMeans2(controlledCounts) +1))
    
    # Compute Z-scores and control for large values
    zScores <- (logControlCounts - rowMeans2(logControlCounts)) / rowSds(logControlCounts)
  }
  else if (!is.matrix(zScores)){
    stop("Provided zScores are not a matrix.")
  }
  if (any(is.infinite(zScores))) {
    # Check for infinite values (should be impossible!)
    stop("Z-score matrix contains infinite values.")
  }
  
  # Transpose zScores if aspect ratio larger than 1
  if (ncol(zScores)/nrow(zScores) > 1){
    zScores <- t(zScores)
  } 
  
  # Perform Singular Value Decomposition (SVD) on the matrix of Z-scores 
  # and extract singular values
  sv <- svd(zScores)$d
  
  # Aspect ratio of the (transposed) count matrix
  numRows <- nrow(zScores)
  numCols <- ncol(zScores)
  beta <- numCols / numRows
  
  # Compute the optimal w(beta)
  coef <- (optimalSVHTCoef(beta) / sqrt(medianMarchenkoPastur(numCols, numRows)))
  
  # compute cutoff
  cutoff <- coef * median(sv)
  
  # compute and return rank
  if (any(sv > cutoff)){
    latentDim <- max(which(sv > cutoff))
  } else {
    warning(paste("Optimal latent space dimension is smaller than 2. Check your count matrix and",
         "verify that all samples have the expected number of counts.",
         "(hist(colSums(counts(ods)))).",
         "For now, the latent space dimension is set to 2.", collapse = "\n"))
    latentDim <- 2} 
  cat("Optimal encoding dimension:", latentDim, "\n")
  return(latentDim)
}

#'
#' Calculate the OHT coefficient 
#' 
#' @noRd
optimalSVHTCoef <- function(beta){ 
  # Calculate lambda(beta)
  sqrt(2 * (beta + 1) + (8 * beta) / (beta + 1 + sqrt(beta^2 + 14 * beta + 1)))
}

#'
#' Calculate the median of the Marchenko-Pastur distribution
#' 
#' Formulas are derived from a robust estimator for the noise 
#' parameter gamma in the model Z~ = Z + gamma * E.
#' Gamma(Z~) = sigma / sqrt(n * mu)) with mu: median of the Marcenko-Pastur distribution.
#' More detailed explanations can be found in Gavish and Donoho 2014.
#' 
#' @noRd
medianMarchenkoPastur <- function(ncol, nrow){ 
  # Compute median of Marchenko-Pastur distribution
  beta <- ncol / nrow
  betaMinus <- (1 - sqrt(beta))^2
  betaPlus <- (1 + sqrt(beta))^2
  
  # Initialize range for upper integral boundary
  lobnd <- copy(betaMinus) 
  hibnd <- copy(betaPlus)
  
  
  while ((hibnd - lobnd) > 0.001){ # Iterate until convergence
    x <- seq(lobnd, hibnd, length.out = 10) # Set range of values for upper integral boundary
    y <- rep(0, length(x))
    for (numCols in 1:length(x)){
      # Approximate integral using Gauss-Kronrod Quadrature
      y[numCols] <- quadgk(dmp, a=betaMinus, b=x[numCols], ndf=nrow, pdim=ncol)
    }  
    
    # Set new boundaries for x that yield the closest results to 0.5
    if (any(y < 0.5)){
      lobnd = max(x[y < 0.5])
    }
    
    if (any(y > 0.5)){
      hibnd = min(x[y > 0.5])
    }
  }
  # If hibnd and lobnd are similar enough, return their mean
  return((hibnd + lobnd) / 2.0)
}
