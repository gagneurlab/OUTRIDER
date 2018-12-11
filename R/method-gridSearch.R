#' 
#' Find the optimal encoding dimension
#' 
#' Finds the optimal encoding dimension for a given data set by running a 
#' grid search based on the provided parameter set.
#'
#' @param ods An OutriderDataSet
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
#'
#' @return The optimal encoding dimension
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
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
#' ods1 <- findEncodingDim(ods, params=encDimSearchParams, 
#'         implementation=implementation)
#' plotEncDimSearch(ods1)
#' 
#' ods2 <- findInjectZscore(ods, zScoreParams=zScoreParams,
#'         encDimParams=encDimSearchParams, implementation=implementation)
#' plotEncDimSearch(ods2)
#' 
#' @rdname findEncodingDim
#' @aliases findEncodingDim, findInjectZscore
#' @export
findEncodingDim <- function(ods, 
                    params=seq(5,min(30,ncol(ods) - 1, nrow(ods) - 1), 2),
                    freq=1E-2, zScore=3, sdlog=log(1.6), lnorm=TRUE, 
                    inj='both', ..., BPPARAM=bpparam()){
    
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
    metadata(ods)[['optimalEncDim']] <- NULL
    metadata(ods)[['optimalEncDim']] <- getBestQ(ods)
    
    counts(ods) <- assay(ods, 'trueCounts')
    validateOutriderDataSet(ods)
    return(ods)
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

#' @rdname findEncodingDim 
#' @export
findInjectZscore <- function(ods, freq=1E-2,
                    zScoreParams=c(seq(1.5, 4, 0.5), 'lnorm'),
                    encDimParams=c(seq(3, 40, 3), seq(45,70, 5), 100, 130, 160),
                    inj='both', ..., BPPARAM=bpparam()){
    encDimParams <- encDimParams[encDimParams < min(dim(ods), nrow(ods))]
    
    FUN <- function(idx, grid, ods, RNDseed, ...){
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
        set.seed(RNDseed)
        findEncodingDim(ods, lnorm=lnorm, zScore=z, params=enc, 
                BPPARAM=SerialParam(), ...)
    }
    
    RNDseed <- .Random.seed
    parGrid <- expand.grid(z=zScoreParams, enc=encDimParams)
    odsres <- bplapply(seq_len(nrow(parGrid)), FUN, 
            ods=ods, freq=freq, inj=inj, RNDseed=RNDseed,
            grid=parGrid, BPPARAM=BPPARAM)
    
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
injectOutliers <- function(ods, freq, zScore, inj, lnorm, sdlog){
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
            counts[idxRow, idxCol] <- art_out
        }else{
            index[idxRow, idxCol] <- 0
            zScore
        }
    }
    # save coruppted counts and index of corruption into ods
    assay(ods, 'counts') <- matrix(as.integer(counts),nrow(ods))
    assay(ods, 'trueCorruptions') <- index
    assay(ods, 'injectedZscore') <- zScore
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
    return(pr$auc.integral)
}

evalAutoCorrection <- function(ods, encoding_dim, BPPARAM=bpparam(), ...){
    
    ods <- OUTRIDER(ods, controlData=TRUE,q=encoding_dim, BPPARAM=BPPARAM, ...)
    eloss <- evalAucPRLoss(ods)
    
    print(paste0('Evaluation loss: ', eloss))
    return(eloss)
}

