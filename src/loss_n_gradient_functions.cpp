// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

const double MIN_EXP_VALUE = -700;

arma::mat minValForExp(arma::mat y){
    arma::uvec idx = find(y < MIN_EXP_VALUE);
    y.elem(idx).fill(MIN_EXP_VALUE);
    return y;
}

// [[Rcpp::export()]]
arma::mat predictMatY(arma::mat x, arma::mat E, arma::mat D, arma::vec b){
    arma::mat y = x * E * D.t();
    y.each_row() += b.t();
    y = minValForExp(y);
    
    return y;
}

// [[Rcpp::export()]]
arma::mat predictMatC(arma::mat x, arma::mat E, arma::mat D, arma::vec b,
                    arma::vec sf){
    arma::mat y = predictMatY(x, E, D, b) ;
    arma::mat c = arma::exp(y);
    c.each_col() %= sf;
    
    return c;
}

arma::vec colMeans(arma::mat X){
    return arma::vectorise(arma::sum(X,0))/X.n_rows;
}

// [[Rcpp::export()]]
double truncLogLiklihoodD(arma::vec par, arma::mat H, arma::vec k, arma::vec sf,
                    arma::vec exclusionMask, double theta, arma::vec thetaC){
    double b, ll, c;
    arma::vec d, y, t1, t2;
    
    arma::uvec idx = find(exclusionMask == 1);
    k = k.elem(idx);
    sf = sf.elem(idx);
    H = H.rows(idx);
    thetaC = thetaC.elem(idx);
    
    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);
    arma::vec thetaVec = theta * thetaC;
    
    y = H * d + b;
    y = minValForExp(y);
    
    t1 = k % (arma::log(sf) + y);
    t2 = (k + thetaVec) % (arma::log(sf) + y + arma::log(1 + thetaVec / (sf % arma::exp(y))));
    
    ll = arma::accu(t1 - t2)/k.n_elem;
    
    return arma::as_scalar(-ll);
}

// [[Rcpp::export()]]
arma::vec gradientD(arma::vec par, arma::mat H, arma::vec k, arma::vec sf,
                    arma::vec exclusionMask, double theta, arma::vec thetaC){
    double b, c;
    arma::vec d, y, yexp, k1, kt, t1, t2, dd, db;
    
    arma::uvec idx = find(exclusionMask == 1);
    k = k.elem(idx);
    sf = sf.elem(idx);
    H = H.rows(idx);
    thetaC = thetaC.elem(idx);
    
    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);
    arma::vec thetaVec = theta * thetaC;
    
    y = H * d + b;
    y = minValForExp(y);
    yexp = arma::exp(y);
    
    t1 = colMeans(k % H.each_col());
    
    kt = (k + thetaVec) / (1 + thetaVec / (sf % yexp));
    t2 = colMeans(kt % H.each_col());
    
    db = arma::vec(1);
    db[0] = arma::accu(kt - k)/k.n_elem;
    dd = t2 - t1;
    
    arma::mat ans = arma::join_cols(db, dd);
    return ans.col(0);
}

// [[Rcpp::export()]]
double truncLogLiklihoodE(arma::vec e, arma::mat D, arma::mat k, arma::vec b,
                        arma::mat x, arma::vec sf, arma::vec theta, 
                        arma::mat exclusionMask, arma::vec thetaC){
    arma::mat E, y, thetaMat, y_plus_log_sf,  t2_1, t2_2, t1, t2, ll;
    
    E = arma::reshape(e, D.n_rows, D.n_cols);
    thetaMat = thetaC * theta.t();
    y = predictMatY(x, E, D, b);
        
    y_plus_log_sf = y.each_col() + arma::log(sf);
    t1 = k % y_plus_log_sf;
    
    t2_2 = arma::exp(y);
    t2_2.each_col() %= sf;
    t2 = (k + thetaMat) % (y_plus_log_sf + arma::log(1 + thetaMat/t2_2));
    
    ll = arma::accu((t1 - t2) % exclusionMask)/arma::accu(exclusionMask);
    return arma::as_scalar(-ll);
}

// [[Rcpp::export()]]
arma::mat gradientE(arma::vec e, arma::mat D, arma::mat k, arma::vec b,
                    arma::mat x, arma::vec sf, arma::vec theta, 
                    arma::mat exclusionMask, arma::vec thetaC){
    arma::mat E, y, thetaMat, t1, kt_2, kt, t3, dE;
    
    E = arma::reshape(e, D.n_rows, D.n_cols);
    thetaMat = thetaC * theta.t();
    y = predictMatY(x, E, D, b);
    
    kt_2 = arma::exp(y);
    kt_2.each_col() %= sf;
    kt = (k + thetaMat) / (1 + thetaMat/kt_2);
    
    t1 = x.t() * ((k  % exclusionMask) * D);
    t3 = x.t() * ((kt % exclusionMask) * D);
    
    dE = (-t1 + t3)/arma::accu(exclusionMask);
    return dE;
}
