debugMyCode <- FALSE

updateD <- function(ods, theta, control, BPPARAM, ...){
    D <- getD(ods)
    b <- getb(ods)
    H <- getx(ods) %*% getE(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    
    fitD <- function(i, D, b, k, H, sf, theta, control){
        pari <- c(b[i], D[i,])
        ki <- k[,i]
        thetai <- theta[i]
        
        fit <- optim(pari, fn=lossDtrunc, gr=gradD, k=ki, H=H, sf=sf, theta=thetai,
                method='L-BFGS', control=control) 
                # lower=rep(-5, ncol(D)+1), upper=rep(5, ncol(D) + 1))
        fit
    }
    
    # TODO check errors: ERROR: ABNORMAL_TERMINATION_IN_LNSRCH
    # This comes from genes where extrem values are present (Z score > 15)
    fitls <- bplapply(1:nrow(ods), fitD, D=D, b=b, k=k, sf=sf, H=H, theta=theta, 
            control=control, BPPARAM=BPPARAM)
    
    # update D and bias terms
    parMat <- sapply(fitls, '[[', 'par')
    print(table(sapply(fitls, '[[', 'message')))
    ods <- setb(ods, parMat[1,])
    ods <- setD(ods, t(parMat)[,-1])
    metadata(ods)[['fits']] <- fitls
    
    return(ods)
}

parametricDFit <- function(){
    ok <- k
    ob <- b
    oD <- D
    oTheta <- theta
    
    i <- 5
    pari <- c(b[i], D[i,])
    ki <- k[,i]
    thetai <- theta[i]
    
    fit <- optim(pari, fn=lossD, gr=gradD, k=ki, H=H, sf=sf, theta=thetai,
                 method='L-BFGS', control=control) 
    # lower=rep(-5, ncol(D)+1), upper=rep(5, ncol(D) + 1))
    fit
    
    coefs <- c(0.1, 1)
    iter <- 0
    while (TRUE) {
        residuals <- disps/(coefs[1] + coefs[2]/means)
        good <- which((residuals > 1e-04) & (residuals < 15))
        suppressWarnings({
            fit <- glm(disps[good] ~ I(1/means[good]), family = Gamma(link = "identity"), 
                       start = coefs)
        })
        oldcoefs <- coefs
        coefs <- coefficients(fit)
        if (!all(coefs > 0)) 
            stop(simpleError("parametric dispersion fit failed"))
        if ((sum(log(coefs/oldcoefs)^2) < 1e-06) & fit$converged) 
            break
        iter <- iter + 1
        if (iter > 10) 
            stop(simpleError("dispersion fit did not converge"))
    }
    names(coefs) <- c("asymptDisp", "extraPois")
    ans <- function(q) coefs[1] + coefs[2]/q
    attr(ans, "coefficients") <- coefs
    ans
    
}

lossD <- function(par, k, H, sf, theta, minMu=0.01){
    b <- par[1]
    d <- par[-1]
    
    y <- H %*% d + b
    yexp <- sf * (minMu + exp(y))
    #yexp <- pmin(1e8, yexp)
    
    ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    
    if(!is.finite(ll) & debugMyCode){
        browser()
    }
    
    return(-ll)
}

lossDtrunc <- function(par, k, H, sf, theta, minMu=0.01){
    b <- par[1]
    d <- par[-1]
    
    y <- H %*% d + b
    #yexp <- sf * (minMu + exp(y))
  
    #ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    # LL = k * log(mu) - (k + theta)*log(mu + theta)  
    t1 <- k * (log(sf) + y + log(1 + minMu/exp(y)))
    t2 <- (k + theta) * (log(sf) + y + log(1 + minMu/exp(y))  + log(1+theta/(sf * (minMu + exp(y))))  )
    ll <- mean(t1 - t2)
    # if(!is.finite(ll) & debugMyCode){
    #   browser()
    # }
  
    return(-ll)
}

lossDtruncDebug <- function(par, k, H, sf, theta, minMu=0.01){
    b <- par[1]
    d <- par[-1]
  
    y <- H %*% d + b
    yexp <- sf * (minMu + exp(y))
  
    #ll <- mean(dnbinom(k, mu=yexp, size=theta, log=TRUE))
    
    ll = mean(k * log(yexp) - (k + theta)*log(yexp + theta)  )
    mean(k*log(yexp))
    mean(k * (log(sf) + y + log(1 + minMu/exp(y))))
    
    # t1 <- k * (log(sf) + y + log(1 + minMu/exp(y)))
    # t2 <- (k + theta) * (log(sf) + y + log(1 + minMu/exp(y))  + log(1+theta/(sf * (minMu + exp(y))))  )
    # ll <- mean(t1 - t2)
    
    # if(!is.finite(ll) & debugMyCode){
    #   browser()
    # }
  
  return(-ll)
}

gradD <- function(par, k, H, sf=1, theta, minMu=0.01){
    b <- par[1]
    d <- par[-1]
    
    y <- c(H %*% d + b)
    yexp <- sf * (minMu + exp(y))
    #yexp <- pmin(1e8, yexp)
    
    k1 <- k * sf * exp(y) / yexp 
    t1 <- colMeans(k1 * H)
    
    kt <- (k + theta) * sf * exp(y) / (yexp + theta)
    t2 <- colMeans(kt * H) 
    
    dd <- t2 - t1
    db <- mean(kt - k1)
    
    
    if(any(!c(is.finite(db), is.finite(dd))) & debugMyCode){
        browser()
    }
    
    return(c(db, dd))
}


debugLossD <- function(){
    samples <- 80
    q<- 2
    
    #Hidden Space is a sample time q matrix. 
    H <- matrix(c(rep(c(-1,1), each=samples/2), 
                  rep(c(-2,2), samples/2)), ncol=2)
    
    H
    D_true <- rnorm(q)
    y_true <- H %*% D_true + 3
    mu_true <- 0.01 + exp(y_true)
    
    k <- rnbinom(length(mu_true), mu = mu_true, size=25)
    
    #library(MASS)
    #glm.nb(k~H)
    
    init<-c(mean(log(k+1)),0,0)
    #b <- 3
    #d <- D_true
    
    
    fit <- optim(init, fn=lossD, gr=gradD, k=k, H=H, s=1, theta=25, method='L-BFGS')
    fit$par
    
    D_true
    lossD(init, k, H, s=1, 25)
    lossD(c(3, D_true), k, H, s=1,  25)
    gradD(c(3, D_true), k, H, s=1,  25)
    gradD(fit$par, k, H, s=1,  25)
    
    numericLossGrad <- function(fn, epsilon, w,...){
        grad <- numeric(length(w))
        for(i in seq_along(w)){
            eps <- integer(length(w))
            eps[i] <- epsilon
            grad[i] <- (fn(w + eps, ...) - fn(w -eps, ...))/(2*epsilon)
        }
        return(grad)
    }
    
    par <- init
    par <- rnorm(11)
    sf<- rnorm(samples, 1, 0.02)
    
    plotNumGrodDiff <- function(){
        numgradloss <- numericLossGrad(lossD, 1E-8, par, k=k, H=H, sf=sf, theta=25)
        gradloss <- gradD(par, k, H, sf=sf, 25)
        
        plot(as.integer((numgradloss - gradloss)*1e7) * (numgradloss - gradloss < 0))
        abline(h=0, col='red')
    }
    plotNumGrodDiff()
    
    sf <- 1
    sf <- sizeFactors(ods)
    par <- pari
    b <- pari[1]
    d <- pari[-1]
    k <- ki
    H <- H
    sf <- sf
    theta <- 0.25
    plot(numericLossGrad(lossD, 1E-8, pari, k=ki, H=H, sf=sf, theta=25),
         gradD(pari, ki, H, sf=sf, 25))
    abline(0,1)
    
    
    which(sapply(fitls, '[[', 'message') == 'ERROR: ABNORMAL_TERMINATION_IN_LNSRCH')
    i <- 437
    i <- 5
    control$trace <- 6
    debugMyCode <- TRUE
    fitD(i=i, D=D, b=b, k=k, sf=sf, H=H, theta=theta, control=control)
    par <- c(0.00940116, 0.0553455, -0.0601702, 0.493793, -0.297911, -0.272072)
    par <- c(mean(log(ki + 1)), D[i,])
    par <- c(b[i], D[i,])
    ki <- k[,i]
    mu <- predictED(getE(ods), getD(ods), getb(ods), getx(ods), sizeFactors(ods))
    cD <- cooksDistance(k, mu, q=5, trim=0.1)
    table(cD > 1)
    sort(round(cD[i,], 2))
    thetai <- theta[i]
    lossD(par, ki, H, sf, theta=thetai)
    gradD(par, ki, H, sf, theta=thetai)
    hist(log10(ki))
    
    sum(ki > 200000)
    sort(ki)
    fit <- optim(c(par), fn=lossD, gr=gradD, k=ki, H=H, sf=sf, theta=thetai,
                 method='L-BFGS', control=control)
    fitD(i=48, D=D, b=b, k=k, sf=sf, H=H, theta=theta, control=control)
    
    lk <- log2(k +1)
    sort((lk - mean(lk))/sd(lk))
    
    lk <- log2(counts(ods) + 1)

}

