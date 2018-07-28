replaceOutliersCooksOutrider <- function(k, mu, q, theta=FALSE, 
                    BPPARAM=bpparam(), ...){
    cooks <- t(cooksDistance(k, mu, q=q))
    ans <- replaceCounts(k, mu, cooks, q=q, ...)
    
    return(list(kReplaced=ans$kReplaced, replacedIndex=ans$index))
}

replaceCounts <- function(k, mu, cooks, q, timMean=0.2, ...){
    
    normFactors <- estimateSizeFactorsForMatrix(t(k))
    if(!missing(mu)){
        normFactors <- mu / exp(colMeans(log(mu)))
    }
    
    # get corrected mean
    ncts <- k/normFactors
    globmean <- apply(ncts, 2, mean, trim=timMean)
    globmean <- matrix(globmean, ncol=ncol(ncts), nrow=nrow(ncts), byrow=TRUE)
    muCorrected <- globmean * normFactors
    
    # get cooks cutof
    cooksCutoff <- qf(0.99, q + 1, nrow(k) - (q + 1))
    idx <- which(cooks > cooksCutoff)
    message(length(idx), ' counts replaced by means.')
    k[idx] <- muCorrected[idx]
    
    return(list(kReplaced=k, index=idx))
}

#'
#' Cooks distance for OUTRIDER
#' Some implementation hints: 
#' \url{http://de.mathworks.com/help/stats/cooks-distance.html}
#' \url{https://en.wikipedia.org/wiki/Cook\%27s_distance}
#' 
#' @noRd
cooksDistance <- function(k, mu, q, timMean=0.2){
    k <- t(k)
    
    # estimate mu if not present
    if(missing(mu)){
        ncts <- t(t(k)/estimateSizeFactorsForMatrix(k))
        timmedMeans <- apply(ncts, 1, mean, trim=trimMean)
        mu <- matrix(trimmedMeans, ncol=ncol(ncts), nrow=nrow(ncts))
    } else {
        mu <- t(mu)
    }
    
    disps <- robustMethodOfMomentsDispOutrider(k, mu)
    
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
robustMethodOfMomentsDispOutrider <- function(cts, mu, minDisp=0.04){
    # get trimmed variance
    v <- DESeq2:::trimmedCellVariance(cts, factor(rep(1, ncol(cts))))
    alpha <- (v - mu)/mu^2
    alpha <- pmax(alpha, minDisp)
    
    return(alpha)
}

