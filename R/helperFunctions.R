
#'
#' This is a helper function to generate a nice color pallet for plotting
#' 
#' @examples
#'   x <- 10
#'   getXColors(x)
#'   
#'   x <- letters[1:10]
#'   getXColors(x)
#' @noRd  
getXColors <- function(x, set="Set2"){
    n <- x
    if(!isScalarNumeric(x)){
        x <- factor(x)
        n <- length(levels(x))
    }
    
    if(length(set) == 1){
        cols <- suppressWarnings(brewer.pal(n=n, name=set))
    } else {
        cols <- set
    }
    
    if(length(cols) > n){
        cols <- cols[seq_len(n)]
    }
    
    if(length(x) > 1){
        ans <- cols[x]
        names(ans) <- as.character(x)
        return(ans)
    }
    
    return(cols)
}


#' This function is used by the plotDispEsts function.
#' 
#' TODO
#' 
#' @noRd
getDispEstsData <- function(ods, mu=NULL){
    if(!'disp' %in% colnames(mcols(ods))){
        stop('Please fit the ods first. ods <- fit(ods)')
    }
    odsMu <- rowMeans(counts(ods, normalized=TRUE))
    if(is.null(mu)){
        mu <- odsMu
    }
    disp <- mcols(ods)$disp
    xidx <- 10^(seq.int(max(-5,log10(min(mu))-1), log10(max(mu))+0.1, 
            length.out = 500))
    
    # fit DESeq2 parametric Disp Fit
    fit <- parametricDispersionFit(mu, 1/disp)
    pred <- fit(xidx)
    return(list(
        mu=mu,
        disp=disp,
        xpred=xidx,
        ypred=pred,
        fit=fit
    ))
}

#'
#' This function is not exported from DESeq2. Therefore we copied it over to
#' here. If DESeq2 will export this function we can use it instead.
#' 
#' TODO
#' 
#' @noRd
parametricDispersionFit <- function (means, disps){
    coefs <- c(0.1, 1)
    iter <- 0
    while (TRUE) {
        residuals <- disps/(coefs[1] + coefs[2]/means)
        good <- which((residuals > 1e-04) & (residuals < 15))
        suppressWarnings({
            fit <- glm(disps[good] ~ I(1/means[good]), 
                    family=Gamma(link="identity"), start=coefs)
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
