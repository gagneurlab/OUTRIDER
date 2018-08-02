

peer <- function(ods, maxItr=1000){
    # check for PEER
    if(!require(peer)){
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
