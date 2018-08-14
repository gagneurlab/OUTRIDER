// https://stackoverflow.com/questions/35923787/fast-large-matrix-multiplication-in-r
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

const double MIN_EXP_VALUE = -500;
    
void printd(SEXP x){
    Rf_PrintValue(x);
}

void printmat(arma::mat x){
    printd(Rcpp::wrap(x.n_cols));
    printd(Rcpp::wrap(x.n_rows));
    
    x = x.head_cols(5);
    x = x.head_rows(5);
    printd(Rcpp::wrap(x));
    return;
}

void pintd(double x){
    printd(Rcpp::wrap(x));
}

arma::mat minValForExp(arma::mat y){
    arma::uvec idx = find(y < MIN_EXP_VALUE);
    y.elem(idx).fill(MIN_EXP_VALUE);
    // y.elem(idx) = arma::vec(idx.n_elem).fill(MIN_EXP_VALUE);
    return y;
}

//we can change this to arma::sum() which sums by col or row.
arma::vec colMeans(arma::mat X){
    arma::vec out = arma::vec(X.n_cols);
    for(int j=0; j < X.n_cols; j++ ) {
        out[j] = arma::accu(X.col(j))/X.n_rows;
    }
    return out;
}

// [[Rcpp::export()]]
double truncLogLiklihoodD(arma::vec par, arma::mat H, arma::vec k, arma::vec sf,
                    arma::vec exclusionMask, double theta, double minMu=0.01,
                    bool verbose=false){
    double b, ll, c;
    arma::vec d, y, ls_y_lmexpy, t1, t2;
    
    arma::uvec idx = find(exclusionMask == 1);
    k = k.elem(idx);
    sf = sf.elem(idx);
    H = H.rows(idx);
    
    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);
    
    y = H * d + b;
    y = minValForExp(y);
    ls_y_lmexpy = arma::log(sf) + y + arma::log(1 + minMu / arma::exp(y));
    
    t1 = k % ls_y_lmexpy;
    t2 = (k + theta) % (ls_y_lmexpy + arma::log(1 + theta / (sf % (minMu + arma::exp(y)))));
    
    ll = arma::accu(t1 - t2)/k.n_elem;
    
    return arma::as_scalar(-ll);
}

// [[Rcpp::export()]]
arma::vec gradientD(arma::vec par, arma::mat H, arma::vec k, arma::vec sf,
                    arma::vec exclusionMask, double theta, double minMu=0.01,
                    bool verbose=false){
    double b, c;
    arma::vec d, y, yexp, k1, kt, t1, t2, dd, db;
    
    arma::uvec idx = find(exclusionMask == 1);
    k = k.elem(idx);
    sf = sf.elem(idx);
    H = H.rows(idx);
    
    b = par.at(0);
    d = par.subvec(1, par.n_elem-1);
    
    if(verbose == true){
        printmat(H);
    }
    
    y = H * d + b;
    y = minValForExp(y);
    yexp = arma::exp(y);
    
    if(verbose == true){
        printmat(y);
    }
    
    k1 = k / (1 + minMu / yexp);
    t1 = colMeans(k1 % H.each_col());
    //printd(Rcpp::wrap(t1));
    
    kt = (k + theta) / (1 + (minMu + theta / sf) / yexp);
    t2 = colMeans(kt % H.each_col());
    //printd(Rcpp::wrap(t2));
    
    db = arma::vec(1);
    db[0] = arma::accu(kt - k1)/k.n_elem;
    dd = t2 - t1;
    
    arma::mat ans = arma::join_cols(db, dd);
    return ans.col(0);
}


// [[Rcpp::export()]]
double truncLogLiklihoodE(arma::vec e, arma::mat D, arma::mat k, arma::vec b,
                        arma::mat x, arma::vec sf, arma::vec theta, 
                        arma::mat exclusionMask, double minMu=0.01){
    arma::mat E, xED, y, ls_y_lmexpy, t1_2, t2_1, t2_2, t1, t2, ll;
    
    E = arma::reshape(e, D.n_rows, D.n_cols);
    xED = x * E * D.t();
    y = xED.each_row() + b.t();
    y = minValForExp(y);
    
    t1_2 = y + arma::log(1 + minMu / arma::exp(y));
    t1_2.each_col() += arma::log(sf);
    t1 = k % t1_2;
    
    t2_1 = k.each_row() + theta.t();
    t2_2 = minMu + arma::exp(y);
    t2_2.each_col() %= sf;
    t2_2 = t1_2 + arma::log(1 + theta.t()/t2_2.each_row());
    t2 = t2_1 % t2_2;
    
    ll = arma::accu((t1 - t2) % exclusionMask)/arma::accu(exclusionMask);
    return arma::as_scalar(-ll);
}


// [[Rcpp::export()]]
arma::mat gradientE(arma::vec e, arma::mat D, arma::mat k, arma::vec b,
                    arma::mat x, arma::vec sf, arma::vec theta, 
                    arma::mat exclusionMask, double minMu=0.01){
    arma::mat E, xED, y, k1, t1, kt_1, kt_2, kt, t3, dE;
    
    E = arma::reshape(e, D.n_rows, D.n_cols);
    xED = x * E * D.t();
    y = xED.each_row() + b.t();
    y = minValForExp(y);
    
    k1 = k / (1 + minMu / arma::exp(y));
    k1 %= exclusionMask;
    t1 = x.t() * (k1 * D);
    
    kt_2 = arma::exp(y);
    kt_2.each_col() %= sf;
    kt_2 = theta.t() / kt_2.each_row();
    kt_2 += 1 + minMu/arma::exp(y);
    kt_1 = k.each_row() + theta.t();
    kt = kt_1 / kt_2;
    kt %= exclusionMask;
    
    t3 = x.t() * (kt * D);
    
    dE = (-t1 + t3)/arma::accu(exclusionMask);
    return dE;
}

// [[Rcpp::export()]]
arma::mat predictY(arma::mat x, arma::mat E, arma::mat D, arma::vec b,
               arma::vec sf, double minMu=0.01){
    
    arma::mat y = x * E * D.t();
    y.each_row() += b.t();
    y = minMu + arma::exp(minValForExp(y));
    y.each_col() %= sf;
    
    return y;
}