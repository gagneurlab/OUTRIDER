if(FALSE){
    library(glmnet)
    devtools::load_all('.')
    ods <- readRDS("../OUTRIDER-analysis-final/Output/data/fitOutrider/ed_NCR_YTC_YRL/Kremer_ODS.RDS")
    ods <- readRDS("../OUTRIDER-analysis-final/Output/data/fitOutrider/ed_NCR_YTC_YRL/SimulationNBinom_ODS.RDS")
    
    ods1 <- autoCorrect(ods, loops=6, implementation = 'ed_NCR_NTC_YRL_YLAS')
    ods2 <- fitAutoencoder(ods2, getBestQ(ods), robust=TRUE, thetaRange=c(0.1,250),
            convergence=1e-5, noRobustLast=FALSE, pValCutoff=0.001, BPPARAM=BPPARAM, 
                   CoxR=FALSE, lasso=TRUE, correctTheta=FALSE, loops=7, initialize = TRUE)
    odsx <- ods2
    odsx <- computePvalues(odsx)
    odsx <- computeZscores(odsx)
    #Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    #0.013     3.148     9.843    34.293    27.023 17620.137 
    par(mfrow=c(2,2))
    plotAberrantPerSample(odsx)
    plotQQ(odsx, global=TRUE, legendPos = NA)
    hist(pValue(odsx))
    par(mfrow=c(1,1))
    plotExpressionRank(odsx, which(rowSums(exclusionMask(odsx) == 0) == 2)[42], normalized = TRUE, basePlot=TRUE)
    plotExpressionRank(odsx, 'CLPP', norm=FALSE)
    i <- 1
    BPPARAM <- MulticoreParam(40, 40, progressbar=TRUE)
    BPPARAM <- SerialParam()
    register(BPPARAM)
    
    tail(sort(aberrant(odsx, by='s')))
    
    hist(lambda(ods2))
    
    abnorm <- mcols(ods)[,'FitDMessage'] != 'CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH'
    lambda(ods)[abnorm]
    plot(sort(counts(ods)[abnorm,][7,]))
    plot(sort(log10(abs(D(ods)[,]))), xlim=c(1, 1000), pch=16, cex=0.5)
    points(sort(log10(abs(D(ods)[abnorm,]))), xlim=c(1, 1000), col='red')
}

fitLassoQ <- function(ods, BPPARAM){
    H <- H(ods)
    k <- counts(ods)
    
    lassoFitLs <- bplapply(1:nrow(ods), singleLassoQ, H=H, k=k, BPPARAM=BPPARAM)
    
    lambda(ods) <- sapply(lassoFitLs, "[", 'lambda.min')
    metadata(ods)[['lassoFits']] <- lassoFitLs
    
    print('Lambda fit')
    print(summary(lambda(ods)))
    
    ods
}

singleLassoQ <- function(i, k, H){
    ki <- k[i,]
    fit <- cv.glmnet(x=H, y=ki, family='poisson')
    return(c(
        lambda.min = fit$lambda.min[1],
        nzero.min  = fit$nzero[fit$lambda == fit$lambda.min][[1]],
        lambda.1se = fit$lambda.1se,
        nzero.1se  = fit$nzero[fit$lambda == fit$lambda.1se][[1]]
    ))
}

