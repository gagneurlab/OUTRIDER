

peer <- function(ods){
    require(peer)
    
    logCts <- log2((counts(ods)+1)/sizeFactors(ods))
    #PEER run.
    model=PEER()
    
    # If no prior information on the
    # magnitude of confounding effects is available, we recommend using 25% of the 
    # number of individuals contained in the study but no more than 100 factors.
    n_unobserved_factors <- min(as.integer(0.25* ncol(ods)),100)
    
    PEER_setNk(model, n_unobserved_factors)
    PEER_setPhenoMean(model, logCts)
    #PEER_setPhenoMean(model, as.matrix(log2(fpkm(ods)+2)))
    PEER_setAdd_mean(model, TRUE)
    
    PEER_update(model)
    
    peerResiduals <- PEER_getResiduals(model)
    peerMean = sizeFactors(ods) * 2^(logCts - peerResiduals)
    
    normalizationFactors(ods) <- pmax(peerMean, 1E-10)
    return(ods)
}