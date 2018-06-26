#' 
#' Find the optimal encoding dimension
#' 
#' Finds the optimal encoding dimension for a given data set by running a 
#' grid search based on the provided params.
#'
#' @param ods An OutriderDataSet
#' @param params Set of possible q values.
#' @param freq Frequency of outlier by default 1E-2
#' @param zScore Injection Z-score, by default 3.
#' @param inj Injection strategy, by default 'both'.
#' @param BPPARAM BPPARAM object by default bpparam().
#'
#' @return The optimal encoding dimension
#' @export
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' findEncodingDim(ods)
findEncodingDim <- function(ods, params=seq(5,30,5), freq=1E-2, zScore=3, 
                    inj='both', BPPARAM=bpparam()){
    
    # compute auto Correction
    ods <- estimateSizeFactors(ods)
    ods <- injectOutliers(ods, freq=freq, zScore=zScore, inj=inj)
    
    ## Check limits of loss:
    # Check min and max of loss by substituting
    # counts and true Counts for the normalization Factors.
    
    # Max overfitting loss:
    normalizationFactors(ods) <- counts(ods) + 1E-9
    overfit <- evalLoss(ods)
    # Best possible loss:
    normalizationFactors(ods) <- assay(ods, 'trueCounts') + 1E-9
    bestloss <- evalLoss(ods)
    # Max underfitting loss (only taking the means of the data):
    normalizationFactors(ods) <- matrix(rep(rowMeans(counts(ods)), 
                                            ncol(ods)),ncol=ncol(ods)) + 1E-9
    underfit <- evalLoss(ods)
    print(paste0('Best possible loss (c=k): ', bestloss,
                '  Max overfitting loss (c=kcorr): ', overfit,
                '  Max underfitting loss (c=colMeans(kcorr): ', underfit))

    #params <- c(5,10,15,20,25,30,50,100)
    eval <- bplapply(params, 
        function(i) evalAutoCorrection(ods, encoding_dim = i),
            BPPARAM=BPPARAM)
    
    eval <- matrix(unlist(eval), ncol = 1, byrow = TRUE)
    colnames(eval) <- c('evalLoss')
    result <- as.data.table(cbind(params,eval))
    
    print(result[which.min(evalLoss),])
    plot(result[,params], result[,evalLoss], log='x', pch=16)
    return(result)
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
    
    #inject counts
    max_out <- max(counts(ods))
    #TODO switch for inj.
    if(inj=='low'){
        index <- -abs(index)
    }
    if(inj=='high'){
        index <- abs(index)
    }
    list_index <- which(index != 0, arr.ind = TRUE)
    # compute size factor normalized counts.
    # don't use it on the ods to not influence the later calculation.
    sf <- estimateSizeFactorsForMatrix(counts(ods))
    normtable <- t(t(counts(ods))/sf)
    counts <- counts(ods)
    
    n_rejected <- 0
    for(i in seq_len(nrow(list_index))){
        cts <- as.numeric(normtable[list_index[i,'row'],])
        fc <- zScore * sd(log2(1 + cts))
        clcount <- index[list_index[i,'row'], list_index[i,'col']]*fc + 
            log2(1 + cts[list_index[i,'col']])
        #multiply size factor again
        art_out <- round(sf[list_index[i,'col']]*2^clcount)
        if(art_out < 100*max_out){
            counts[list_index[i,'row'], list_index[i,'col']] <- art_out
        }else{
            #remove super large outliers
            index[list_index[i,'row'], list_index[i,'col']] <- 0 
            n_rejected <- n_rejected + 1
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

evalAutoCorrection <- function(ods, encoding_dim=20, theta=25){
    
    ods <- autoCorrect(ods, q=encoding_dim)
    eloss <- evalLoss(ods, theta)
    
    print(paste0('Evaluation loss: ', eloss))
    return(eloss)
}

