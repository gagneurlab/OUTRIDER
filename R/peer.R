

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


runSva <- function(ods){
    ods <- readRDS("./../scared-analysis/Output/data/GTEx_not_sun_exposed_OutriderDONE.RDS")
    ods <- makeExampleOutriderDataSet(dataset="Kremer")[1:140, 1:140]
    k <- counts(ods)
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
    
    
    k <- counts(ods)
    
    ####
    # prepare data
    ###
    dim(k)
    
    # random group assignment
    group = sample(c(TRUE, FALSE), ncol(k), replace = TRUE)
    sum(group)/length(group)
    
    # create model
    dat0 = as.matrix(log10(k+1))
    mod1 = model.matrix(~group)
    mod0 = cbind(mod1[,1])
    
    # estimate q
    n.sv <- num.sv(dat0, mod0)
    n.sv
    
    # run sva
    svseq = svaseq(dat0,mod1,mod0,n.sv=15)
    
    svseq
    plot(svseq,pch=19,col="blue")
    
    library(bladderbatch)
    data(bladderdata)
    pheno = pData(bladderEset)
    
    contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
    fitContrasts = contrasts.fit(fit,contrast.matrix)
    

    return(ods)
}
