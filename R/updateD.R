#'
#' Update D function
#' 
#' @noRd
updateD <- function(ods, minMu, control, BPPARAM){
    D <- D(ods)
    b <- b(ods)
    H <- x(ods) %*% E(ods)
    k <- t(counts(ods))
    sf <- sizeFactors(ods)
    mask <- exclusionMask(ods)
    theta <- theta(ods)
    
    fitls <- bplapply(1:nrow(ods), singleDFit, D=D, b=b, k=k, sf=sf, H=H, 
            theta=theta, mask=mask, minMu=minMu, control=control, 
            BPPARAM=BPPARAM)
    
    # update D and bias terms
    parMat <- sapply(fitls, '[[', 'par')
    print(table(sapply(fitls, '[[', 'message')))
    mcols(ods)[,'NumConvergedD'] <- mcols(ods)[,'NumConvergedD'] + grepl(
        "CONVERGENCE: REL_REDUCTION_OF_F .. FACTR.EPSMCH", 
        sapply(fitls, '[[', 'message'))
    b(ods) <- parMat[1,]
    D(ods) <- t(parMat)[,-1]
    
    metadata(ods)[['Dfits']] <- fitls
    
    return(ods)
}


singleDFit <- function(i, D, b, k, theta, mask, ...){
    pari <- c(b[i], D[i,])
    ki <- k[,i]
    thetai <- theta[i]
    maski <- mask[i,]
    
    fit <- optim(pari, fn=truncLogLiklihoodD, gr=gradientD, 
            k=ki, theta=thetai, exclusionMask=maski, ...,
            lower=-100, upper=100, method='L-BFGS')
    return(fit)
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
    yexp <- sf * exp(y)
  
    #ll = mean(k * log(yexp) - (k + theta)*log(yexp + theta))

    t1 <- k * (log(sf) + y + log(1 + minMu/exp(y)))
    t2 <- (k + theta) * (log(sf) + y + log(1 + minMu/exp(y))  + log(1+theta/(sf * (minMu + exp(y))))  )
    ll <- mean(t1 - t2)

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
    
    #k1 <- k * sf * exp(y) / yexp 
    k1 <- k / (1 + minMu/exp(y) )
    t1 <- colMeans(k1 * H)
    
    #kt <- (k + theta) * sf * exp(y) / (yexp + theta)
    kt <- (k + theta) / ( 1 + (minMu + theta/sf)/exp(y) )
    
    t2 <- colMeans(kt * H) 
    
    dd <- t2 - t1
    db <- mean(kt - k1)
    
    
    if(any(!c(is.finite(db), is.finite(dd))) & debugMyCode){
        browser()
    }
    
    return(c(db, dd))
}
