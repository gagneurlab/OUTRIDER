// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

arma::mat MatMultAtBC(arma::mat A, arma::mat B, arma::mat C){
    arma::mat At = A.t();
    arma::mat D = At * B * C;
    return(D);
}

arma::mat predict(arma::mat x, arma::mat W, arma:: mat b){
    arma::mat D = x * W * W.t();
    arma::mat y = (D + b);
    return(y);
}

arma::mat computeLoss(arma::mat k, arma::mat y, arma::mat s, arma::mat theta){
    
    arma::mat sexpy = s % arma::exp(y);
    arma::mat m1=-arma::sum(k % (arma::log(s) + y),1);
    arma::mat m2= arma::sum((k + theta) % arma::log(sexpy + theta),1);
    arma::mat m3= arma::sum(m1+m2, 0);
    return(m3);
}


arma::mat computeKT(arma::mat k, arma::mat x, arma::mat W, arma::mat b, arma::mat s, double theta){
    arma::mat thetaMat(size(k));
    thetaMat.fill(theta);
    arma::mat y = predict(x, W, b);
    arma::mat sexpy = s % arma::exp(y);
    arma::mat kt = (k + thetaMat)%sexpy/(sexpy+thetaMat);
    return(kt);
}


//TODO try mat A(5, 6);
// A.fill(123.0); to set theta.
// [[Rcpp::export]]
SEXP truncLogLiklihood(arma::mat k, arma::mat x, arma::mat W, arma::mat b, arma::mat s, double theta){
    arma::mat thetaMat(size(k));
    thetaMat.fill(theta);
    arma::mat y = predict(x, W, b);
    
    arma::mat Loss = computeLoss(k, y, s, thetaMat);
    return Rcpp::wrap(Loss);
}

// [[Rcpp::export]]
SEXP gradLogLiklihood(arma::mat k, arma::mat x, arma::mat W, arma::mat b, arma::mat s, double theta){
    
    arma::mat kt = computeKT(k, x, W, b, s, theta);
    
    arma::mat t1 = MatMultAtBC(x, k, W);
    arma::mat t2 = MatMultAtBC(k, x, W);
    
    arma::mat t3 = MatMultAtBC(x, kt, W);
    arma::mat t4 = MatMultAtBC(kt, x, W);
    
    arma::mat dw = (-t1 - t2 + t3 + t4);
    arma::mat db = arma::sum(kt-k,0).t();
    arma::mat grad = join_rows(dw, db);
    return Rcpp::wrap(grad);
}


