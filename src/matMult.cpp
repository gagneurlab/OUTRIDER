// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
    
// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;
    
    return Rcpp::wrap(C);
}

// function for 3 matrices.
// [[Rcpp::export]]
SEXP armaMatMultAtBC(arma::mat A, arma::mat B, arma::mat C){
    arma::mat At = A.t();
    arma::mat D = At * B * C;
    
    return Rcpp::wrap(D);
}
// [[Rcpp::export]]
SEXP armaMatMultABBt(arma::mat A, arma::mat B){
    arma::mat Bt = B.t();
    arma::mat D = A * B * Bt;
    
    return Rcpp::wrap(D);
}
