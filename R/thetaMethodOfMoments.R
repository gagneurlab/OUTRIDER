#' 
#' Method of moments for theta calculations
#' 
#' Have a look for more details at 
#' \code{\link[DESeq2]{robustMethodOfMomentsDisp}}
#'
#
#' This function is adapted from the DESeq2:::trimmedCellVariance
#' to only one factor. This reduces some computation
#'   
#' @noRd
robustMethodOfMomentsOfTheta <- function(cts, maxTheta, minTheta, minMu=0.01){
    
    mue <- apply(cts, 1, mean, trim=1/8)
    mue <- pmax(mue, minMu)
    se <- (cts - mue)^2
    see <- apply(se, 1, mean, trim=1/8)
    ve <- 1.51*see
    
    theta <- mue^2/(ve-mue)
    
    theta[theta < 0] <- maxTheta
    theta <- pmax(minTheta, pmin(maxTheta, theta))
    return(theta)
}
