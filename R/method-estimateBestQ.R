#' 
#' Find the optimal encoding dimension
#' 
#' Finds the optimal encoding dimension by either Optimal Hard Thresholding 
#' or injecting artificial splicing outlier ratios while maximizing the 
#' precision-recall curve.
#'

#'
#' @return The optimal encoding dimension
#'

#' @param ods An OutriderDataSet object
#' @param zScoresOHT A z-score matrix
#' @param useOHT If \code{TRUE} (default), Optimal Hard Thresholding is 
#'              used to estimate the optimal encoding dimension.
#' @param params,encDimParams Set of possible q values.
#' @param freq Frequency of outlier, defaults to 1E-2
#' @param zScore,zScoreParams Set of possible injection Z-score, defaults to 3.
#' @param sdlog Standard deviation of the sitribution on the log scale.
#' @param lnorm If TRUE, the default, Z-scores are drawn from a log normal 
#'             distribution with a mean of \code{log(zScore)} in log-scale.
#' @param inj Injection strategy, by default 'both'.
#' @param ... Further arguments passed on to the \code{controlForConfounders}
#'             function.
#' @param BPPARAM BPPARAM object by default bpparam().
#' @import RMTstat
#' @import pracma
#' @return The OutriderDataSet object with the optimal encoding dimension saved in the metadata
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' 
#' # run OHT (default)
#' estimateBestQ(ods)
#' 
#' # run hyperparameter optimization (grid-search)
#' encDimSearchParams <- c(5, 8, 10, 12, 15)
#' zScoreParams <- c(2, 3, 5, 'lnorm')
#' implementation <- 'autoencoder'
#' register(MulticoreParam(4))
#' \dontshow{
#'     ods <- ods[1:12,1:12]
#'     encDimSearchParams <- c(2)
#'     zScoreParams <- c('lnorm')
#'     register(SerialParam())
#'     implementation <- 'pca'
#' }
#' ods1 <- estimateBestQ(ods, useOHT=FALSE, params=encDimSearchParams, 
#'         implementation=implementation)
#' plotEncDimSearch(ods1)
#' 
#' ods2 <- findInjectZscore(ods, zScoreParams=zScoreParams,
#'         encDimParams=encDimSearchParams, implementation=implementation)
#' plotEncDimSearch(ods2) 
#' 
#' @rdname estimateBestQ
#' @aliases estimateBestQ, findInjectZscore
#' @export
estimateBestQ <- function(ods=NULL, zScoresOHT=NULL, useOHT=TRUE, 
                          params=seq(2, min(100, ncol(ods) - 1, nrow(ods) - 1), 2),
                          freq=1E-2, zScore=3, sdlog=log(1.6), lnorm=TRUE,
                          inj='both', ..., BPPARAM=bpparam()){
  # Optimal Hard Thresholding (default)
  if (isTRUE(useOHT)){
    if (is.null(ods) & is.null(zScoresOHT)){
      stop("Please provide an OutriderDataSet or a z-score matrix.")
    }
    else if (!is.null(ods)){
      OUTRIDER:::checkOutriderDataSet(ods)
      if (!is.null(zScoresOHT)){
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
      zScoresOHT <- (logControlCounts - rowMeans2(logControlCounts)) / rowSds(logControlCounts)
    }
    else if (!is.matrix(zScoresOHT)){
      stop("Provided zScoresOHT are not a matrix.")
    }
    if (any(is.infinite(zScoresOHT))) {
      # Check for infinite values (should be impossible!)
      stop("Z-score matrix contains infinite values.")
    }
    
    # Transpose zScoresOHT if aspect ratio larger than 1
    if (ncol(zScoresOHT)/nrow(zScoresOHT) > 1){
      zScoresOHT <- t(zScoresOHT)
    } 
    
    # Perform Singular Value Decomposition (SVD) on the matrix of Z-scores 
    # and extract singular values
    sv <- svd(zScoresOHT)$d
    
    # Aspect ratio of the (transposed) count matrix
    numRows <- nrow(zScoresOHT)
    numCols <- ncol(zScoresOHT)
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
                    "verify that all samples have the expected number of counts",
                    "(hist(colSums(counts(ods)))).",
                    "For now, the latent space dimension is set to 2.", collapse = "\n"))
      latentDim <- 2} 
    cat("Optimal encoding dimension:", latentDim, "\n")
    
    if (!is.null(ods)){
      data <- data.table(singular_values=sv)
      data <- data[, q:=.I]
      data <- data[-.N] # remove last entry (close to zero, impedes visualization)
      data <- data[, oht:=FALSE]
      data <- data[q==latentDim, oht:=TRUE] # set optimal latent dimension
      data <- data[order(-oht)]
      
      metadata(ods)[['encDimTable']] <- NULL
      metadata(ods)[['encDimTable']] <- data
      metadata(ods)[['optimalEncDim']] <- latentDim
      
      validateOutriderDataSet(ods)
      return(ods)
    } else{
      return(latentDim)
    }
    
  } else{ 
    # Hyperparamter optimization (grid-search)
    
    # compute auto Correction
    ods <- estimateSizeFactors(ods)
    ods <- injectOutliers(ods, freq=freq, zScore=zScore, inj=inj, lnorm=lnorm,
                          sdlog=sdlog)
    
    # TODO Please check if this is still correct
    # ods <- getLossBounderyCases(ods)
    
    eval <- bplapply(X=params, ..., BPPARAM=BPPARAM, 
                     FUN=function(i, ..., evalAucPRLoss=NA){
                       evalAutoCorrection(ods, encoding_dim=i, BPPARAM=SerialParam(), ...)}
    )
    
    metadata(ods)[['encDimTable']] <- data.table(
      encodingDimension= params,
      evaluationLoss= unlist(eval), 
      evalMethod='aucPR')
    metadata(ods)[['optimalEncDim']] <- getBestQ(ods)
    
    counts(ods) <- assay(ods, 'trueCounts')
    validateOutriderDataSet(ods)
    return(ods)
  }
}

getLossBounderyCases <- function(ods){
    
    ## Check limits of loss:
    # Check min and max of loss by substituting
    # counts and true Counts for the normalization Factors.
    
    # Max overfitting loss:
    normalizationFactors(ods) <- pmax(counts(ods), 1E-9)
    overfit <- evalLoss(ods)
    
    # Best possible loss:
    normalizationFactors(ods) <- pmax(assay(ods, 'trueCounts'), 1E-9)
    bestloss <- evalLoss(ods)
    
    # Max underfitting loss (only taking the means of the data):
    normalizationFactors(ods) <- matrix(
                rep(pmax(rowMeans(counts(ods)), 1E-9), ncol(ods)),
            ncol=ncol(ods))
    underfit <- evalLoss(ods)
    
    metadata(ods)[['boundaryCases']] = c(
            'Best possible loss (c=k)' = bestloss,
            'Max overfitting loss (c=kcorr)' = overfit,
            'Max underfitting loss (c=colMeans(kcorr)' = underfit)
    
    return(ods)
}

#' @rdname estimateBestQ
#' @export
findInjectZscore <- function(ods, freq=1E-2,
                    zScoreParams=c(seq(1.5, 4, 0.5), 'lnorm'),
                    encDimParams=c(seq(3, 40, 3), seq(45,70, 5), 100, 130, 160),
                    inj='both', ..., BPPARAM=bpparam()){
    encDimParams <- encDimParams[encDimParams < min(dim(ods), nrow(ods))]
    
    FUN <- function(idx, grid, ods, RNGseed, ...){
        z <- grid[idx,"z"]
        enc <- grid[idx,"enc"]
        lnorm <- FALSE
        if(z == 'lnorm'){
            z <- 3
            lnorm=TRUE
        } else {
            z <- as.integer(z)
        }
        message(date(), ": run Z-score: ", z, " and enc: ", enc)
        set.seed(RNGseed)
        estimateBestQ(ods, useOHT=FALSE, lnorm=lnorm, zScore=z, params=enc, 
                BPPARAM=SerialParam(), ...)
    }
    
    RNGseed <- .Random.seed
    parGrid <- expand.grid(z=zScoreParams, enc=encDimParams)
    odsres <- bplapply(seq_row(parGrid), FUN, ods=ods, freq=freq,
            inj=inj, grid=parGrid, RNGseed=RNGseed, BPPARAM=BPPARAM)
    
    res <- rbindlist(lapply(seq_along(odsres), function(i){
            metadata(odsres[[i]])$encDimTable[,.(
                    encodingDimension, zScore=parGrid[i,"z"], evaluationLoss)]
    }))
    res <- rbindlist(lapply(zScoreParams, function(z){
        tmpdt <- res[zScore==z]
        tmpdt[,opt:=getBestQDT(tmpdt)]
        tmpdt
    }))
    metadata(ods)[['encDimTable']] <- res
    return(ods)
}

## Helper functions called within estimateBestQ ##


#' injectOutliers
#'
#' @param ods OutriderDataSet
#' @param freq frequency of injected counts.
#' @param zScore injection Z-score.
#' @param inj injection strategy.
#'
#' @return and OutriderDataSet with artificially corrupted counts.
#' @noRd
injectOutliers <- function(ods, freq, zScore, inj, lnorm, sdlog){
    # copy true counts to be able to acces them in the loss later
    assay(ods, 'trueCounts', withDimnames=FALSE) <- counts(ods)
    
    # generate index of injected counts
    size <- prod(dim(ods))
    index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2), 
            replace = TRUE)
    index <- matrix(index, nrow = nrow(ods))
    switch(inj,
        low = { index <- -abs(index) },
        high = { index <- abs(index) }
    )
    
    tmpzScore <- matrix(0, ncol=ncol(ods), nrow=nrow(ods))
    if(isTRUE(lnorm)){
        tmpzScore[index!=0] <- rlnorm(sum(index!=0), log(zScore), sdlog=sdlog)
    } else {
        tmpzScore[index!=0] <- zScore
    }
    zScore <- tmpzScore
    
    #inject counts
    max_out <- min(10*max(counts(ods), na.rm=TRUE), .Machine$integer.max)
    
    # compute size factor normalized counts.
    # don't use it on the ods to not influence the later calculation.
    sf <- estimateSizeFactorsForMatrix(counts(ods))
    normtable <- t(t(counts(ods))/sf)
    counts <- counts(ods)
    
    list_index <- which(index != 0, arr.ind = TRUE)
    for(i in seq_len(nrow(list_index))){
        idxCol <- list_index[i,'col']
        idxRow <- list_index[i,'row']
        
        cts <- as.numeric(normtable[idxRow,])
        fc <- zScore[idxRow, idxCol] * sd(log2(1 + cts))
        clcount <- index[idxRow, idxCol] * fc + log2(1 + cts[idxCol])
        
        #multiply size factor again
        art_out <- round(sf[idxCol] * 2^clcount)
        
        # only insert outliers if they are different from before 
        # and not too large
        if(art_out < max_out & counts[idxRow, idxCol] != art_out){

          # do not introduce outlier such that all samples have 0 counts for that gene.
          # if all other samples (excluding idxCol) are 0 and art_out is 0 do not inject that outlier.
          if(all(counts[idxRow,-idxCol] == 0) & art_out == 0){
              index[idxRow, idxCol] <- 0
          }else{
              counts[idxRow, idxCol] <- art_out
          }
        }else{
            index[idxRow, idxCol] <- 0
            zScore
        }
    }
    # save coruppted counts and index of corruption into ods
    assay(ods, 'counts', withDimnames=FALSE) <- matrix(as.integer(counts), 
            nrow=nrow(ods))
    assay(ods, 'trueCorruptions', withDimnames=FALSE) <- index
    assay(ods, 'injectedZscore', withDimnames=FALSE) <- zScore
    return(ods)
}


#' 
#' evaluate loss = 1/n * Sum(log-likelihood(k | c))
#' 
#' @noRd
evalLoss <- function(ods, theta=25){
    N_corrupted <- sum(abs(assay(ods, 'trueCorruptions')))
    kTrue <- assay(ods, 'trueCounts')[assay(ods, 'trueCorruptions')!=0]
    c <- assay(ods, 'normalizationFactors')[assay(ods, 'trueCorruptions')!=0]
    eval <- sum(-dnbinom(kTrue, mu=c , size=theta, log=TRUE))/N_corrupted
    return(eval)
}

evalAucPRLoss <- function(ods){
    scores <- -as.vector(assay(ods, 'pValue'))
    labels <- as.vector(assay(ods, 'trueCorruptions') != 0) + 0
    
    if(any(is.na(scores))){
        warning(sum(is.na(scores)), " P-values where NAs.")
        scores[is.na(scores)] <- min(scores, na.rm=TRUE)-1
    }
    pr <- pr.curve(scores, weights.class0=labels)
    return(max(0, pr$auc.integral, na.rm=TRUE))
}

evalAutoCorrection <- function(ods, encoding_dim, BPPARAM=bpparam(), ...){
    
    ods <- OUTRIDER(ods, controlData=TRUE,q=encoding_dim, BPPARAM=BPPARAM, ...)
    eloss <- evalAucPRLoss(ods)
    
    print(paste0('Evaluation loss: ', eloss,' for q=',encoding_dim))
    return(eloss)
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

