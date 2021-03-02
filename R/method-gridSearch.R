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
findEncodingDim <- function(ods, prepro_options=getDefaultPreproParams(ods),
                    params=getParamsToTest(ods),
                    freq=1E-2, zScore=3, sdlog=log(1.6), lnorm=TRUE,
                    inj='both', ..., BPPARAM=bpparam()){
    
    # compute auto Correction
    ods <- preprocess(ods, prepro_options)
    ods <- injectOutliers.OUTRIDER2(ods, freq=freq, zScore=zScore, inj=inj, 
                                    lnorm=lnorm, sdlog=sdlog)
    # no preprocessing after injecting outliers as outiers are already injected
    # into preprocessed values
    prepro_options$prepro_func <- NULL 
    
    dot_args <- list(...)
    if("usePython" %in% names(dot_args) & isTRUE(dot_args[["usePython"]]) ){
        eval <- bplapply(X=params, ..., BPPARAM=SerialParam(), 
                        FUN=function(i, ..., evalAucPRLoss=NA){
                            evalAutoCorrection(ods, encoding_dim=i, 
                                            prepro_options=prepro_options, 
                                            BPPARAM=BPPARAM, ...)}
        )
    } else{
        eval <- bplapply(X=params, ..., BPPARAM=BPPARAM, 
            FUN=function(i, ..., evalAucPRLoss=NA){
                evalAutoCorrection(ods, encoding_dim=i, 
                                    prepro_options=prepro_options, 
                                    BPPARAM=SerialParam(), ...)}
        )
    }
    
    metadata(ods)[['encDimTable']] <- data.table(
            encodingDimension= params,
            evaluationLoss= unlist(eval), 
            evalMethod='aucPR')
    metadata(ods)[['optimalEncDim']] <- NULL
    metadata(ods)[['optimalEncDim']] <- getBestQ(ods)
    
    preprocessed(ods) <- assay(ods, 'trueObservations') # TODO check 
    validateOutrider2DataSet(ods)
    return(ods)
}

#' @rdname findEncodingDim 
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
        findEncodingDim(ods, lnorm=lnorm, zScore=z, params=enc, 
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
    assay(ods, 'trueObservations', withDimnames=FALSE) <- counts(ods)
    
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
    assay(ods, 'counts', withDimnames=FALSE) <- matrix(as.integer(counts), 
            nrow=nrow(ods))
    assay(ods, 'trueCorruptions', withDimnames=FALSE) <- index
    assay(ods, 'injectedZscore', withDimnames=FALSE) <- zScore
    return(ods)
}

#' injectOutliers
#'
#' @param ods OutriderDataSet
#' @param freq frequency of injected counts.
#' @param zScore injection Z-score.
#' @param inj injection strategy.
#'
#' @return and OutriderDataSet with artificially corrupted counts.
#' @noRd
injectOutliers.OUTRIDER2 <- function(ods, freq, zScore, inj, lnorm, sdlog){
    checkOutrider2DataSet(ods)
    checkPreprocessing(ods)
    
    # copy true values to be able to acces them in the loss later
    assay(ods, 'trueObservations', withDimnames=FALSE) <- observed(ods)
    
    # generate index of injected values
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
    
    #inject values
    if(is.integer(preprocessed(ods))){
        max_out <- min(10*max(observed(ods), na.rm=TRUE), .Machine$integer.max)
    } else{
        max_out <- min(10*max(observed(ods), na.rm=TRUE), .Machine$double.xmax)
    }
    
    sf <- sizeFactors(ods) # equal to 1 if not requested
    normtable   <- transformed(ods) # already sf normalized (if requested)
    revtransFUN <- metadata(ods)$prepro_options$reverse_data_trans
    if(is.null(revtransFUN)){
        revtransFUN <- identity
    } else{
        revtransFUN <- match.fun(revtransFUN)
    }
    
    values <- preprocessed(ods)
    
    list_index <- which(index != 0, arr.ind = TRUE)
    for(i in seq_len(nrow(list_index))){
        idxCol <- list_index[i,'col']
        idxRow <- list_index[i,'row']
        
        x_trans <- as.numeric(normtable[idxRow,])
        fc <- zScore[idxRow, idxCol] * sd(x_trans, na.rm=TRUE)
        x_inj <- index[idxRow, idxCol] * fc + x_trans[idxCol]
        
        #multiply size factor again
        art_out <- getReverseTransformed(x_trans=x_inj, sf=sf[idxCol], 
                                            revtransFUN=revtransFUN)
        if(profile(ods) == "outrider"){
            art_out <- round(art_out)
        }
        
        # only insert outliers if they are different from before 
        # and not too large and ensure values are >= 0 if profile is outrider 
        # or protrider
        if(!is.na(values[idxRow, idxCol]) & 
                values[idxRow, idxCol] != art_out &
                art_out < max_out & 
                art_out < 1000 * max(values[idxRow,], na.rm=TRUE) & 
                !(art_out < 0 && (profile(ods) != "other")) ){
            values[idxRow, idxCol] <- art_out
        }else{
            index[idxRow, idxCol] <- 0
            zScore[idxRow, idxCol] <- 0
        }
    }
    
    if(profile(ods) == "outrider"){
        values <- as.integer(values)
    }
    
    # save coruppted counts and index of corruption into ods
    observed(ods, withDimnames=FALSE) <- matrix(values, nrow=nrow(ods))
    assay(ods, 'trueCorruptions', withDimnames=FALSE) <- index
    assay(ods, 'injectedZscore', withDimnames=FALSE) <- zScore
    metadata(ods)$prepro_options$prepro_func <- NULL
    return(ods)
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
    
    ods <- OUTRIDER(ods, controlData=TRUE, q=encoding_dim,
                    BPPARAM=BPPARAM, ...)
    # ods <- OUTRIDER(ods, q=encoding_dim, BPPARAM=BPPARAM, ...)
    eloss <- evalAucPRLoss(ods)
    
    print(paste0('Evaluation loss: ', eloss,' for q=',encoding_dim))
    return(eloss)
}

getParamsToTest <- function(ods, maxTestedDimensionProportion=3){
    # q parameters to be tested (from DROP)
    a <- 5 
    b <- min(ncol(ods), nrow(ods)) / maxTestedDimensionProportion   # N/3
    
    maxSteps <- 15
    if(maxTestedDimensionProportion < 4){
        maxSteps <- 20
    }
    
    Nsteps <- min(maxSteps, b)   # Do at most 20 steps or N/3
    # Do unique in case 2 were repeated
    pars_q <- unique(round(exp(seq(log(a),log(b),length.out = Nsteps))))
    
    return(pars_q)
}
