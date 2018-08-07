replaceOutliersCooksOutrider <- function(k, mu, q, theta=FALSE, 
                    BPPARAM=bpparam(), ...){
    cooks <- t(cooksDistance(k, mu, q=q))
    ans <- replaceCounts(k, mu, cooks, q=q, ...)
    
    return(list(kReplaced=ans$kReplaced, replacedIndex=ans$index))
}

replaceCounts <- function(k, mu, cooks, q, trim=0.2, ...){
    
    normFactors <- estimateSizeFactorsForMatrix(t(k))
    if(!missing(mu)){
        normFactors <- mu / exp(colMeans(log(mu)))
    }
    
    # get corrected mean
    ncts <- k/normFactors
    globmean <- apply(ncts, 2, mean, trim=trim)
    globmean <- matrix(globmean, ncol=ncol(ncts), nrow=nrow(ncts), byrow=TRUE)
    muCorrected <- globmean * normFactors
    
    # get cooks cutof
    cooksCutoff <- qf(0.99, q + 1, nrow(k) - (q + 1))
    idx <- which(cooks > cooksCutoff)
    message(length(idx), ' counts replaced by means.')
    k[idx] <- round(muCorrected[idx])
    
    return(list(kReplaced=k, index=idx))
}

#'
#' Cooks distance for OUTRIDER
#' Some implementation hints: 
#' \url{http://de.mathworks.com/help/stats/cooks-distance.html}
#' \url{https://en.wikipedia.org/wiki/Cook\%27s_distance}
#' 
#' @noRd
cooksDistance <- function(k, mu, q, trim=0.2){
    k <- t(k)
    
    # estimate mu if not present
    if(missing(mu)){
        ncts <- t(t(k)/estimateSizeFactorsForMatrix(k))
        trimmedMeans <- apply(ncts, 1, mean, trim=trim)
        mu <- matrix(trimmedMeans, ncol=ncol(ncts), nrow=nrow(ncts))
    } else {
        mu <- t(mu)
    }
    
    disps <- robustMethodOfMomentsOFDispersion(k, mu)
    
    # calculate H (based on DESeq2)
    w <- (mu^-1 + disps)^-1
    xtwx <- rowSums(w)
    H <- w * xtwx^-1;
    
    V <- mu + disps * mu^2
    PearsonResSq <- (k - mu)^2/V
    
    cooksD <- (PearsonResSq/(q + 1)) * H/(1 - H)^2
    return(cooksD)
}


#' 
#' Method of moments for dispersion calculations
#' 
#' Have a look for more details at 
#' \code{\link[DESeq2]{robustMethodOfMomentsDisp}}
#' 
#' @noRd
robustMethodOfMomentsOFDispersion <- function(cts, mu, maxDisp=100){
    v <- trimmedVariance(cts)
    theta <- mu^2/(v-mu)
    theta <- pmin(theta, maxDisp)
    return(theta)
}

#'
#' This function is adapted from the DESeq2:::trimmedCellVariance
#' to only one factor. This reduces some computation
#' 
#' @noRd
trimmedVariance <- function(cts){
    mue <- apply(cts, 1, mean, trim=1/8)
    se <- (cts - mue)^2
    see <- apply(se, 1, mean, trim=1/8)
    ve <- 1.51*see
    return(ve)
}


estimateThetaFromCounts <- function(k, mu, mask, BPPARAM=bpparam()){
    ods <- OutriderDataSet(countData=k)
    ods <- estimateSizeFactors(ods)
    if(missing(mu)){
        mu <- matrix(rowMeans(counts(ods, normalized=TRUE)), 
                nrow=nrow(ods), ncol=ncol(ods))
    }
    normalizationFactors(ods) <- mu
    ods <- fit(ods, excludeMask=mask, BPPARAM=BPPARAM)
    return(ods)
}

findOutlierNBfit<- function(k, mu, pValCutoff=0.001){
    k <- t(k)
    mu <- t(mu)
    
    ods <- OutriderDataSet(countData=k)
    normalizationFactors(ods) <- mu
    ods <- fit(ods)
    ods <- computePvalues(ods)
    mask <- assay(ods, 'pValue') < pValCutoff
    message(sum(mask),' outliers excluded from fit.')
    ods <- fit(ods, excludeMask=mask)
    theta <- mcols(ods)[['disp']]
    
    mask <- t(mask)
    ans <- list(mask=mask)
    ans[['theta']]<- theta
    return(ans)
}
