

replaceOutliersCooksOutrider <- function(k, mu, q, theta=FALSE, 
                                 BPPARAM=bpparam(), ...){
    cooks <- t(cooksDistance(k, mu, q=q))
    
    kReplaced <- replaceCounts(k,mu,cooks, q=q, ...)
}

replaceCounts <- function(k, mu, cooks, q, ...){
    cooksCutoff <- qf(0.99, q, nrow(k) - q)
    
    if(missing(mu)){
        s <- estimateSizeFactorsForMatrix(t(k))
        ncts <- k/s
        globmean <- apply(ncts, 2, mean, trim=0.2)
        globmean <- t(matrix(globmean, ncol=ncol(ncts), nrow=nrow(ncts), byrow=TRUE))
        muCorrected <- globmean * s
    }else{
        normFactors <- mu / exp(colMeans(log(mu)))
        ncts <- k/normFactors
        globmean <- apply(ncts, 2, mean, trim=0.2)
        globmean <- matrix(globmean, ncol=ncol(ncts), nrow=nrow(ncts), byrow=TRUE)
        muCorrected <- globmean * normFactors
    }
    
    idx <- which(cooks > cooksCutoff)
    message(length(idx), ' counts replaced by means.')
    k[idx] <- muCorrected[idx]
    return(k)
}


cooksDistance <- function(k, mu, w, q){
    k <- t(k)
    #w <- w_fit
    # estimate mu if not present
    if(missing(mu)){
        ncts <- t(t(k)/estimateSizeFactorsForMatrix(k))
        mu <- matrix(trimmedMean(ncts), ncol=ncol(ncts), nrow=nrow(ncts))
    }else{
        mu <- t(mu)
    }
    
    dispersionso <- robustMethodOfMomentsDispOutrider(k, mu)
    
    # calculate H
    # no W from PCA/AutoEncoder 
    if(missing(w)){
        w <- (mu^-1 + dispersionso)^-1
        xtwx <- rowSums(w)
        H <- w * xtwx^-1;
        # W form AutoEncoder/PCA
    } 
    # else {
    #     W <- getWeights(w, nrow(k))
    #     b <- getBias(w, nrow(k))
    #     
    #     # y = h*W.t() + b
    #     # h = x*W
    #     # c = s * exp(y+xbar)
    #     X <- t(k) %*% W %*% t(W)
    #     X <- t(W)
    #     nxbar <- rowMeans(log((k+1)/sf))
    #     nx <- log((k+1)/sf) - nxbar
    #     y = t(nx) %*% W %*% t(W) + b
    #     nc <- sf * exp(y+nxbar)
    #     H <- leverageCalc(t(y))
    # }
    Vo <- mu + dispersionso * mu^2
    PearsonResSqo <- (k - mu)^2/Vo
    
    cookso <- (PearsonResSqo/(q + 1)) * H/(1 - H)^2
    return(cookso)

    # # D = e^2 / (s^2*p) * H/(1-H)^2)
    # # s^2 = (n - p)^-1* e^T
    # e <- k - mu
    # 
    # firstNom <- e^2
    # firstDenom <- as.numeric(solve(ncol(k) - (q+1))) * e^2
    # secTerm <- H/(1 - H)^2
    # 
    # D = (firstNom / firstDenom) %*% secTerm
    # 
    # return(D)
    # 
    # http://de.mathworks.com/help/stats/cooks-distance.html
    # https://en.wikipedia.org/wiki/Cook%27s_distance
}

levCalcR <- function(X){
    # H = X ( X^T * X)^-1 * X^T
    X %*% solve(t(X) * X) * t(X)
}

robustMethodOfMomentsDispOutrider <- function(cts, mu, minDisp=0.04){
    
    # get trimmed variance
    v <- DESeq2:::trimmedCellVariance(cts, factor(rep(1, ncol(cts))))
    
    alpha <- (v - mu)/mu^2
    alpha <- pmax(alpha, minDisp)
    
    return(alpha)
}

trimmedMean <- function(cts, cut=0.1){
    ncut <- floor(ncol(cts)*cut/2)
    vcut <- rep(c(FALSE, TRUE, FALSE), c(ncut, ncol(cts) - ncut*2, ncut))
    
    trmean <- apply(cts, 1, vcut=vcut, function(x, vcut){
        mean(sort(x)[vcut])})
    return(trmean)
}
