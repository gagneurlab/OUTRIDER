#' Low Rank Simulation
#'
#' This function simulates a data set as a low-rank signal corrupted by 
#' Gaussian noise. This function is copied from 
#' https://github.com/julierennes/denoiseR/blob/master/R/LRsim.R as it was 
#' removed from CRAN due to lack of maintenance.
#'  
#' @param n integer, number of rows
#' @param p integer, number of columns
#' @param k integer, rank of the signal
#' @param SNR numeric, signal to noise ratio
#' @return X the simulated data 
#' @return mu the true signal 
#' @return sigma the standard deviation of the noise added to the signal
#' @details A data set of size n*p and of rank k is simulated. More precisely, 
#'   it is simulated as follows: A SVD is performed on a n*p matrix generated
#'   from a standard multivariate normal distribution. Then, the signal is 
#'   computed using the first k singular vectors and singular values 
#'   U_q D_q V_q'. The signal is scaled in such a way that the variance of each
#'   column is 1 and then a Gaussian noise with variance sigma^2 is added.
#'   The SNR is calculated as 1/ sigma sqrt(np). 
#' @examples  
#' Xsim= LRsim(100, 30, 2, 2)
#' 

.LRsim <- function(n, p, k, SNR){
  noisesd <- 1/(SNR*sqrt(n*p)) 
  if(k == 0){
    MU <- matrix(0, n, p)  
  } else {
    SIGNAL <- replicate(p, rnorm(n, 0, 1))
    SIGNAL <- scale(SIGNAL, scale = FALSE)
    svdSIGNAL <- svd(SIGNAL) 
    MU <- (
      svdSIGNAL$u[, 1:k,drop=F] 
      %*% diag(svdSIGNAL$d[1:k], k, k) 
      %*% t(svdSIGNAL$v[, 1:k, drop=F]) 
      / sqrt(sum(svdSIGNAL$d[1:k]^2))
    )
  }
  X <- MU + noisesd*replicate(p, rnorm(n, 0, 1))
  return(list(X = X, mu = MU, sigma = noisesd))
}
