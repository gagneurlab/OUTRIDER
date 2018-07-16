#' 
#' Find the optimal encoding dimension
#' 
#' Finds the optimal encoding dimension for a given data set by running a 
#' grid search based on the provided parameter set.
#'
#' @param ods An OutriderDataSet
#' @param params Set of possible q values.
#' @param freq Frequency of outlier, defaults to 1E-2
#' @param zScore Injection Z-score, defaults to 3.
#' @param inj Injection strategy, by default 'both'.
#' @param BPPARAM BPPARAM object by default bpparam().
#'
#' @return The optimal encoding dimension
#'
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet(180, 60)
#' ods <- findEncodingDim(ods, params=2:9)
#' 
#' # plot the results of the dimension search
#' metadata(ods)$encDimTable[plot(encodingDimension, evaluationLoss)]
#' 
#' @export
findEncodingDim <- function(ods, params=seq(5,min(30,ncol(ods), nrow(ods)), 2),
                    freq=1E-2, zScore=3, inj='both', ..., BPPARAM=bpparam()){
    
    # compute auto Correction
    ods <- estimateSizeFactors(ods)
    ods <- injectOutliers(ods, freq=freq, zScore=zScore, inj=inj)
    
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

    eval <- bplapply(X=params, ..., BPPARAM=BPPARAM, FUN=function(i, ...){
        evalAutoCorrection(ods, encoding_dim=i, BPPARAM=SerialParam(), ...)})
    
    metadata(ods)[['encDimTable']] <- data.table(
            'encodingDimension' = params,
            'evaluationLoss' = unlist(eval))
    
    metadata(ods)[['optimalEncDim']] <- getBestQ(ods)
    
    counts(ods) <- assay(ods, 'trueCounts')
    validateOutriderDataSet(ods)
    return(ods)
}

## Helper functions called within SearchEncDim ##


#' injectOutliers
#'
#' @param ods OutriderDataSet
#' @param freq frequency of injected counts.
#' @param zScore injection Z-score.
#' @param inj injection strategy.
#'
#' @return and OutriderDataSet with artificially corrupted counts.
#' @noRd
injectOutliers <- function(ods, freq, zScore, inj){
    # copy true counts to be able to acces them in the loss later
    assays(ods)[['trueCounts']] <- counts(ods)
    
    # generate index of injected counts
    size <- prod(dim(ods))
    index <- sample(c(0,1,-1), size, prob = c(1 - freq, freq/2, freq/2), 
            replace = TRUE)
    index <- matrix(index, nrow = nrow(ods))
    switch(inj,
        low = { index <- -abs(index) },
        high = { index <- abs(index) }
    )
    
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
        fc <- zScore * sd(log2(1 + cts))
        clcount <- index[idxRow, idxCol] * fc + log2(1 + cts[idxCol])
        
        #multiply size factor again
        art_out <- round(sf[idxCol] * 2^clcount)
        
        # only insert outliers if they are different from before 
        # and not too large
        if(art_out < max_out & counts[idxRow, idxCol] != art_out){
            counts[idxRow, idxCol] <- art_out
        }else{
            index[idxRow, idxCol] <- 0 
        }
    }
    # save coruppted counts and index of corruption into ods
    assay(ods, 'counts') <- matrix(as.integer(counts),nrow(ods))
    assays(ods)[['trueCorruptions']] <- index 
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

evalAutoCorrection <- function(ods, encoding_dim=20, theta=25, ...){
    
    ods <- autoCorrect(ods, q=encoding_dim, ...)
    eloss <- evalLoss(ods, theta)
    
    print(paste0('Evaluation loss: ', eloss))
    return(eloss)
}

