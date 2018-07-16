

peer <- function(ods){
    require(peer)
    
    logCts <- log2(counts(ods)+1)
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
    peerMean = 2^(logCts - peerResiduals)-1
    
    normalizationFactors(ods) <- pmax(peerMean, 1E-10)
    return(ods)
}


runSva <- function(ods){
    ods <- readRDS("./../scared-analysis/Output/data/GTEx_not_sun_exposed_OutriderDONE.RDS")
    ods <- makeExampleOutriderDataSet(dataset="Kremer")[1:140, 1:140]
    k <- counts(ods)
    keep = apply(k, 1, function(x) length(x[x>5])>=2)
    table(keep)
    filtk <- k[keep,]
    genes = rownames(filtk)[grep("^ENS", rownames(filtered))]
    controls <- sample(c(TRUE, FALSE), ncol(k), replace = TRUE)
    table(controls)
    
    library(zebrafishRNASeq)
    data(zfGenes)
    table(grepl("^ERCC", rownames(zfGenes)))
    table(grepl("^ENS", rownames(zfGenes)))
    
    controls = grepl("^ERCC", rownames(filtered))
    
    group = sample(c(TRUE, FALSE), ncol(k), replace = TRUE)
    sum(group)/length(group)
    dat0 = as.matrix(filtk)
    mod1 = model.matrix(~group)
    mod0 = cbind(mod1[,1])
    
    # estimate q
    n.sv <- num.sv(dat0, mod1)
    
    # run sva
    svseq = svaseq(dat0,mod1,mod0,n.sv=15)
    
    svseq
    plot(svseq,pch=19,col="blue")
    
    return(ods)
}