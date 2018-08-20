#' 
#' Method of moments for theta calculations
#' 
#' Have a look for more details at 
#' \code{\link[DESeq2]{robustMethodOfMomentsDisp}}
#' 
#' @noRd
robustMethodOfMomentsOfTheta <- function(cts, mu, maxTheta=250, minTheta=0.1){
    v <- trimmedVariance(cts)
    theta <- mu^2/(v-mu)
    
    theta <- pmax(minTheta, pmin(maxTheta, rowMaxs(theta)))
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
