#'
#' PEER implementation
#'
#' @noRd
peer <- function(ods, maxItr=1000){
    
    # check for PEER
    if(!requireNamespace("peer")){
        stop("Please install the 'peer' package from GitHub to use this ",
                "functionality.")
    }
    checkSizeFactors(ods)
    
    # log counts
    logCts <- log2(t(t(counts(ods)+1)/sizeFactors(ods)))
    
    # default and recommendation by PEER: min(0.25*n, 100)
    n_unobserved_factors <- min(as.integer(0.25* ncol(ods)),100)
    
    # prepare PEER model
    model <- PEER()
    PEER_setNmax_iterations(model, maxItr)
    PEER_setNk(model, n_unobserved_factors)
    PEER_setPhenoMean(model, logCts)
    PEER_setAdd_mean(model, TRUE)
    
    # run fullpeer pipeline
    PEER_update(model)
    
    # extract PEER data
    peerResiduals <- PEER_getResiduals(model)
    peerMean = t(t(2^(logCts - peerResiduals)) * sizeFactors(ods))
    
    # save model in object
    normalizationFactors(ods) <- pmax(peerMean, 1E-8)
    metadata(ods)[["PEER_model"]] <- list(
            alpha     = PEER_getAlpha(model),
            residuals = PEER_getResiduals(model),
            W         = PEER_getW(model))
    
    return(ods)
}

