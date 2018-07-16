#' 
#' Autoencoder function to correct for confounders.
#' 
#' This is the wrapper function for the autoencoder implementation. 
#' It can be used to call the standard R implementation or the experimental
#' Python implementation.
#'
#' @param ods An OutriderDataSet object
#' @param q The encoding dimensions
#' @param theta The dispersion parameter
#' @param implementation "R", the default, will use the R implementation or 
#'             "python" to use the python/tensorflow experimental implementation
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam} instance
#'             to be used for parallel computing.
#' @param ... passed on to the autoencoder implementing method. In the case of 
#'             the R implementation it is passed to the optim function. 
#' 
#' @return An ods object including the control factors 
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' 
#' plotCountCorHeatmap(ods, normalized=FALSE)
#' plotCountCorHeatmap(ods, normalized=TRUE)
#' 
#' @export
autoCorrect <- function(ods, q, theta=25, 
                    implementation=c("R", "python", "PEER", "robustR",
                            "robustRM1","robustRTheta"),
                    BPPARAM=bpparam(), ...){
    
    # error checking
    checkOutriderDataSet(ods)
    checkCountRequirements(ods)
    checkSizeFactors(ods)
    
    if(!missing(q)){
        if(!is.numeric(q) && q > 0){
            stop("Please provide an integer greater then 0 for q.")
        }
        if(q >= nrow(ods)){
            stop("Please use a q smaller than the number of features.")
        }
        if(q >= ncol(ods)){
            stop("Please use a q smaller than the number of samples.")
        }
    } else {
        q <- getBestQ(ods)
        if(is.na(q)){
            q <- 5
        }
    }
    
    # pass on to the correct implementation
    switch(match.arg(implementation),
        R = {
            impl <- "standard R"
            ans <- autoCorrectR(ods, q, theta, ...)
        },
        PEER = {
            impl <- "PEER"
            ans <- peer(ods)
        },
        robustR = {
            impl <- "robust R"
            ans <- autoCorrectRCooksIter2(ods, q, theta, BPPARAM=BPPARAM, ...)
        },
        robustRM1 = {
            impl <- 'robustRM1'
            ans <- autoCorrectRCooksIter(ods, q, theta, ...)
        },
        robustRTheta = {
            impl <- 'robustRTheta'
            ans <- autoCorrectRCooksIterTheta(ods, q, theta, ...)
        },
        python = {
            impl <- "TensorFlow"
            ans <- autoCorrectPython(ods, ...)
        },
        stop("Requested autoCorrect implementation is unknown.")
    )
    
    message(date(), ": Used the ", impl, " implementation for autoCorrect.")
    return(ans)
}

#' 
#' Autoencoder function to correct for confounders.
#'
#' @param ods An uormalized OUTRIDER data set
#' @param q the encoding dimension used.
#' @param theta value used in the likelihood (default=25).
#' 
#' @noRd
autoCorrectR <- function(ods, q, theta=25, control=list(), ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }

    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
                loss(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    fit <- optim(w_guess, loss, gr = lossGrad, k=k, x=x, s=s, xbar=xbar, 
            theta=theta, method="L-BFGS-B", control=control, ...)
    #Check that fit converged
    if(fit$convergence!=0){
        warning(paste0("Fit didn't converge with warning: ", fit$message))
    }
    
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
                loss(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}


#' 
#' Extracting the latent space
#' 
#' Extracts the latent space from the OutriderDataSet object 
#' determined by the autoencoder.
#'
#' @param ods An OutriderDataSet
#'
#' @return A matrix containing the latent space determined by the autoencoder.
#'
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- autoCorrect(ods)
#' computeLatentSpace(ods)[,1:6]
#' 
#' @export
computeLatentSpace <- function(ods){
    stopifnot(is(ods, 'OutriderDataSet'))
    if(!'weights' %in% names(metadata(ods))){
        stop('No weights are stored in this OutriderDataSet object. ',
                'Please compute weights before extracting the latent space.')
    }
    if(any(metadata(ods)[['dim']]!=dim(ods))){
        stop('The OutriderDataSet dimension changed and does not match with ',
                'the existing weights. Please recompute the weights. ',
                'Computation not possible on this data.')
    }
    
    # get data
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    w <- metadata(ods)[['weights']]
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    l <- t(x%*%W) 
    
    if(ncol(l)!=ncol(ods)){
        stop('Dimensions do not match.')
    }
    return(l)
}



## loss and gradient of loss function accesed by the functions above. ##

#' loss function
#'
#' @param w weight matrix 
#' @param k counts
#' @param x log-centered counts
#' @param s size factors
#' @param xbar offset 
#' @param theta value used in the likelihood (default=25).
#'
#' @return Returns the negative log likelihood of the negative binomial 
#' @noRd
loss <- function(w, k, x, s, xbar, theta){
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    
    b <- b + xbar
    truncLogLiklihood(k, x, W, b, s, theta)
}


#' gradient of loss function
#'
#' @param w weight matrix 
#' @param k counts
#' @param x log-centered counts
#' @param s size factors
#' @param xbar offset 
#' @param theta value used in the likelihood (default=25).
#'
#' @return returns the gradient of the loss function.
#' @noRd
lossGrad <- function(w, k, x, s, xbar, theta){
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    
    b <- b + xbar
    gradLogLiklihood(k, x, W, b, s, theta)
}



#'
#' Predict controlled Counts.
#' 
#' @param w weight matrix 
#' @param k counts
#' @param s size factors
#' @param xbar offset 
#'
#' @return Returns the predicted corrections (predicted means).
#' @noRd
predictC <- function(w, k, s, xbar){
    x <-  t(t(log((1+k)/s)) - xbar)
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    
    y <- t(t(x%*%W %*% t(W)) +b + xbar)
    s*exp(y)
}

#'
#' Get the bias term from the weights vector
#' 
#' @noRd
getBias <- function(w, nr){
    w[seq_len(nr)+(length(w)/nr-1)*nr]
}

#'
#' Get the weights matrix from the weight vector without the bias term
#' 
#' @noRd
getWeights <- function(w, nr){
    matrix(w, nrow=nr)[,seq_len(length(w)/nr-1)]
}





autoCorrectR2 <- function(ods, q, theta=25, control=list(), ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    # compute log of per gene centered counts 
    x0 <- log((1+k)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    
    theta <- matrix(
            1/(0.1 + 1/colMeans(k)),
            ncol = ncol(k),
            nrow = nrow(k),
            byrow = TRUE)
    
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
               loss2(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    fit <- optim(w_guess, loss2, gr = lossGrad2, k=k, x=x, s=s, xbar=xbar, 
                 theta=theta, method="L-BFGS-B", control=control, ...)
    #Check that fit converged
    if(fit$convergence!=0){
        warning(paste0("Fit didn't converge with warning: ", fit$message))
    }
    
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               loss2(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}


loss2 <- function(w, k, x, s, xbar, theta){
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    ## encoding 
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    
    y <- t(t(x%*%W %*% t(W)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)
    
    ## log likelihood 
    ll <- dnbinom(k, mu= y_exp, size=theta, log=TRUE)
    - mean( ll )
}




lossGrad2 <- function(w, k, x, s, xbar, theta){
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    
    # dW:
    t1 <- t(x) %*% (k %*% W)
    #t1 <- armaMatMultAtBC(x, k, W)
    t2 <- t(k) %*% (x %*% W)
    #t2 <- armaMatMultAtBC(k, x, W)
    y <- t(t(x%*%W %*% t(W)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)
    kt <- (k + theta)*y_exp/(y_exp+theta)
    t3 <- t(x) %*% (kt %*% W)
    #t3 <- armaMatMultAtBC(x, kt, W)
    t4 <- t(kt) %*% (x %*% W)
    #t4 <- armaMatMultAtBC(kt, x, W)
    dw <- (-t1 - t2 + t3 + t4)/prod(dim(k))
    
    #db:
    db <- colSums(kt-k)/prod(dim(k))
    
    return(c(dw, db))
}

loss3 <- function(w, k, x, s, xbar, theta){
    ## log, size factored, and centered counts 
    #x <-  t(t(log((1+k)/s)) - xbar)
    ## encoding 
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    
    y <- t(t(x%*%W %*% t(W)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y) -1
    
    ## log likelihood 
    ll <- dnbinom(k, mu= pmax(1E-10,y_exp), size=theta, log=TRUE)
    - mean( ll )
}




lossGrad3 <- function(w, k, x, s, xbar, theta){
    W <- matrix(w, nrow=ncol(k))
    b <- W[,ncol(W)]
    W <- W[,1:ncol(W)-1]
    
    # dW:
    y <- t(t(x%*%W %*% t(W)) + xbar + b)
    #y <- t(t(armaMatMultABBt(x, W)) + xbar + b)
    y_exp <- s*exp(y)
    
    k1 <- k*y_exp /pmax(y_exp-1,1E-10)
    k1[which(y_exp-1<1E-10)] <- 0 
    t1 <- t(x) %*% (k1 %*% W)
    #t1 <- armaMatMultAtBC(x, k, W)
    t2 <- t(k1) %*% (x %*% W)
    #t2 <- armaMatMultAtBC(k, x, W)
    
    #kt <- (k + theta)*(y_exp)/(y_exp - 1 + theta)
    kt <- (k + theta)*(y_exp)/(pmax(y_exp-1,1E-10) + theta)
    kt[which(y_exp-1<1E-10)] <- 0 
    t3 <- t(x) %*% (kt %*% W)
    #t3 <- armaMatMultAtBC(x, kt, W)
    t4 <- t(kt) %*% (x %*% W)
    #t4 <- armaMatMultAtBC(kt, x, W)
    dw <- (-t1 - t2 + t3 + t4)/prod(dim(k))
    
    #db:
    db <- colSums(kt-k1)/prod(dim(k))
    
    return(c(dw, db))
}

predictC3 <- function(w, k, s, xbar){
    x <-  t(t(log((1+k)/s)) - xbar)
    b <- getBias(w, ncol(k))
    W <- getWeights(w, ncol(k))
    
    y <- t(t(x%*%W %*% t(W)) +b + xbar)
    s*exp(y)
}

replaceOutliersCooks <- function(k, mu, theta=FALSE, thetaOUTRIDER=TRUE, 
                    BPPARAM=bpparam()){
    k <- t(k)
    dds <- DESeqDataSetFromMatrix(countData=k, design=~1, 
            colData=DataFrame(seq_len(ncol(k))))
    dds <- estimateSizeFactors(dds)
    if(!missing(mu)){
        mu <- t(mu)
        normFactors <- mu/ exp(rowMeans(log(mu)))
        normalizationFactors(dds)  <- normFactors
    }
    dds <- DESeq2:::DESeqParallel(dds, test="Wald", fitType="mean",
            quiet=FALSE, modelMatrix = NULL, useT=FALSE, minmu=0.1, 
            betaPrior=FALSE, BPPARAM=BPPARAM)
    dds <- replaceOutliers(dds)

    if(theta==FALSE){
        return(t(counts(dds)))
    }
    if(isTRUE(thetaOUTRIDER)){
        dds <- OutriderDataSet(dds)
        dds <- fit(dds)    
        return(list(t(counts(dds)), mcols(dds)[['disp']]))
    }else{
        return(list(t(counts(dds)), 1/dispersions(dds)))
    }
    return(t(counts(dds)))
}

autoCorrectRCooksIter <- function(ods, q, theta=25, control=list(), 
                    BPPARAM=BPPARAM, ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    
    k_no <-replaceOutliersCooks(k, BPPARAM=BPPARAM)
    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
               loss(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    
    w_fit <- w_guess
    for(i in 1:10){
        
        k_no <-replaceOutliersCooks(k,predictC(w_fit, k, s, xbar))
        x0 <- log((1+k_no)/s)
        x <- t(t(x0) - xbar)
        
        control$maxit <- 10    
        fit <- optim(w_fit, loss3, gr = lossGrad3, k=k_no, x=x, s=s, xbar=xbar, 
                     theta=theta, method="L-BFGS-B", control=control, ...)
        
        w_fit <- fit$par
        message('Iteration ', i, ' loss: ', loss(w_fit, k, x, s, xbar, theta))
        
        #Check that fit converged
        if(fit$convergence!=0){
            warning(paste0("Fit didn't converge with warning: ", fit$message))
        }
        
    }
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               loss(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC3(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}

autoCorrectRCooksIter2 <- function(ods, q, theta=25, control=list(), 
                    BPPARAM=bpparam(), ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    
    
    k_no <-replaceOutliersCooks(k, BPPARAM=BPPARAM)
    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
               loss(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    
    w_fit <- w_guess
    for(i in 1:10){
        
        k_no <-replaceOutliersCooks(k,predictC(w_fit, k, s, xbar), 
                BPPARAM=BPPARAM)
        x0 <- log((1+k_no)/s)
        x <- t(t(x0) - xbar)
        
        control$maxit <- 10    
        fit <- optim(w_fit, loss, gr = lossGrad, k=k_no, x=x, s=s, xbar=xbar, 
                 theta=theta, method="L-BFGS-B", control=control, ...)
        
        w_fit <- fit$par
        message('Iteration ', i, ' loss: ', loss(w_fit, k, x, s, xbar, theta))
        
        #Check that fit converged
        if(fit$convergence!=0){
            warning(paste0("Fit didn't converge with warning: ", fit$message))
        }
        
    }
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               loss(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    validObject(ods)
    return(ods)
}   


autoCorrectRCooksIterTheta <- function(ods, q, theta=25, control=list(), 
                    delta, downweight, BPPARAM=BPPARAM, ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)

    dispfit <-replaceOutliersCooks(k, theta=TRUE, BPPARAM=BPPARAM)
    k_no <- dispfit[[1]]
    theta <- dispfit[[2]]

    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
               loss2(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    
    w_fit <- w_guess

    for(i in 1:20){
        
        dispfit <-replaceOutliersCooks(k, predictC(w_fit, k = k, s=s, xbar=xbar), theta=TRUE, BPPARAM=BPPARAM)
        k_no <- dispfit[[1]]
        if(i > 2){
            theta <- rowMeans(cbind(theta , dispfit[[2]]))
        }
        message('min: ',min(theta),
                ' 25%: ',quantile(theta, p=0.25),
                ' 50%: ',quantile(theta, p=0.5),
                ' 75%: ',quantile(theta, p=0.75),
                ' max: ',max(theta))
        x0 <- log((1+k_no)/s)
        x <- t(t(x0) - xbar)
        
        control$maxit <- 10    
        fit <- optim(w_fit, loss2, gr = lossGrad2, k=k_no, x=x, s=s, xbar=xbar, 
                     theta=theta, method="L-BFGS-B", control=control, ...)
        
        w_fit <- fit$par
        message('Iteration ', i, ' loss: ', loss2(w_fit, k, x, s, xbar, theta))
        
        #Check that fit converged
        if(fit$convergence!=0){
            warning(paste0("Fit didn't converge with warning: ", fit$message))
        }
        
    }
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               loss2(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    
    mcols(ods)[['thetaAutoenc']] <- theta
    
    validObject(ods)
    return(ods)
}   


autoCorrectRCooksIterThetaM1 <- function(ods, q, theta=25, control=list(), ...){
    
    if(!'factr' %in% names(control)){
        control$factr <- 1E9
    }
    
    k <- t(counts(ods, normalized=FALSE))
    s <- sizeFactors(ods)
    theta <- rep(theta, ncol(k))
    
    dispfit <-replaceOutliersCooks(k, theta=TRUE)
    k_no <- dispfit[[1]]
    #theta <- dispfit[[2]]
    # compute log of per gene centered counts 
    x0 <- log((1+k_no)/s)
    xbar <- colMeans(x0)
    x <- t(t(x0) - xbar)
    
    
    # initialize W using PCA and bias as zeros.
    pca <- pca(x, nPcs = q) 
    pc  <- loadings(pca)
    w_guess <- c(as.vector(pc), numeric(ncol(k)))
    # check initial loss
    print(
        paste0('Initial PCA loss: ',
               loss3(w_guess, k, x, s, xbar, theta))
    )
    
    # optimize log likelihood
    t <- Sys.time()
    
    w_fit <- w_guess
    for(i in 1:10){
        
        dispfit <-replaceOutliersCooks(k, predictC(w_fit, k = k, s=s, xbar=xbar), theta=TRUE)
        k_no <- dispfit[[1]]
        if(i > 2){
            theta <- rowMeans(cbind(theta, dispfit[[2]]))
        }
        message('min: ',min(theta),
                ' 25%: ',quantile(theta, p=0.25),
                ' 50%: ',quantile(theta, p=0.5),
                ' 75%: ',quantile(theta, p=0.75),
                ' max: ',max(theta))

        x0 <- log((1+k_no)/s)
        x <- t(t(x0) - xbar)
        
        control$maxit <- 10    
        fit <- optim(w_fit, loss3, gr = lossGrad3, k=k_no, x=x, s=s, xbar=xbar, 
                     theta=theta, method="L-BFGS-B", control=control, ...)
        
        w_fit <- fit$par
        message('Iteration ', i, ' loss: ', loss3(w_fit, k, x, s, xbar, theta))
        
        #Check that fit converged
        if(fit$convergence!=0){
            warning(paste0("Fit didn't converge with warning: ", fit$message))
        }
        
    }
    w_fit <- fit$par
    print(Sys.time() - t)
    print(
        paste0('nb-PCA loss: ',
               loss2(w_fit,k, x, s,xbar, theta))
    )
    
    correctionFactors <- t(predictC3(w_fit, k, s, xbar))
    stopifnot(identical(dim(counts(ods)), dim(correctionFactors)))
    
    # add it to the object
    normalizationFactors(ods, replace=TRUE) <- correctionFactors
    metadata(ods)[['weights']] <- w_fit
    metadata(ods)[['dim']] <- dim(ods)
    
    mcols(ods)[['thetaAutoenc']] <- theta
    
    validObject(ods)
    return(ods)
}   
