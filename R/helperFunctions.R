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
        if (!all(coefs > 0)){
            warning("Parametric dispersion fit failed.", 
                    " Using last working coefficients:", 
                    paste0(round(oldcoefs, 3), sep=", "))
            coefs <- oldcoefs
            break
        }
        if ((sum(log(coefs/oldcoefs)^2) < 1e-06) & fit$converged) 
            break
        iter <- iter + 1
        if (iter > 100) {
            warning("Dispersion fit did not converge after 100 ",
                    "iterations. We stopped here.")
            break
        }
    }
    names(coefs) <- c("asymptDisp", "extraPois")
    ans <- function(q) coefs[1] + coefs[2]/q
    attr(ans, "coefficients") <- coefs
    ans
}

#'
#' Get the gene name or index 
#' 
#' @noRd
getGeneIndex <- function(geneIdx, ods){
    if(is.null(geneIdx)){
        stop('Please provide a geneID')
    }
    if(is.logical(geneIdx)){
        geneIdx <- which(geneIdx)
    }
    if(is.numeric(geneIdx)){
        if(!(is.numeric(geneIdx) && max(geneIdx) <= nrow(ods))){
            stop('Gene index is out of bounds:', paste(geneIdx, collapse=", "))
        }
        if(!is.null(rownames(ods))){
            geneIdx <- rownames(ods)[geneIdx]
        }
    }
    if(is.character(geneIdx) & any(!geneIdx %in% rownames(ods))){
        stop("Gene ID is not in the data set.")
    }
    return(geneIdx)
}

getBestQ <- function(ods){
    if('optimalEncDim' %in% names(metadata(ods))){
        return(metadata(ods)[['optimalEncDim']])
    }
    if('encDimTable' %in% names(metadata(ods))){
        return(metadata(ods)[['encDimTable']][
            which.min(evaluationLoss),encodingDimension])
    }
    # warning('Please find the optimal encoding dimension by running. ')
    return(NA_integer_)
}

checkCountRequirements <- function(ods, test=FALSE){
    checkOutriderDataSet(ods)
    if(ncol(ods) == 0){
        stop("Please provide at least one sample.")
    }
    if(nrow(ods) == 0){
        stop("Please provide at least one gene.")
    }
    filterGenes <- rowSums(counts(ods)) == 0
    if(any(filterGenes) & isFALSE(test)){
        stop("There are genes without any read. Please filter first ",
                "the data with: ods <- filterExpression(ods)")
    }
    filterGenes <- filterGenes | rowSums(counts(ods)) < ncol(ods)/100
    if(any(filterGenes) & isFALSE(test)){
        stop("The model requires for each gene at least 1 read ", 
                "every 100 sample. Please filte first the data with: ",
                "ods <- filterExpression(ods)")
    }
    return(invisible(filterGenes))
}

checkOutriderDataSet <- function(ods){
    if(!is(ods, 'OutriderDataSet')){
        stop('Please provide an OutriderDataSet')
    }
    return(invisible(TRUE))
}

checkSizeFactors <- function(ods, funName=sys.calls()[[1]]){
    checkOutriderDataSet(ods)
    if(is.null(sizeFactors(ods))){
        stop("Please calculate the size factors before calling the '", funName,
                "' function. Please do: ods <- estimateSizeFactors(ods)")
    }
    return(invisible(TRUE))
}

checkFullAnalysis <- function(ods, funName=sys.calls()[[1]]){
    checkSizeFactors(ods)
    if(!'padjust' %in% assayNames(ods)){
        stop("Please calculate the P-values before calling the '", funName,
             "' function. Please do: ods <- computePvalues(ods)")
    }
    if(!'zScore' %in% assayNames(ods)){
        stop("Please calculate the Z-scores before calling the '", funName,
             "' function. Please do: ods <- computeZscores(ods)")
    }
    return(invisible(TRUE))
}

