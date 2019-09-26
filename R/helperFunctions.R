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
    if(all(is.na(set))){
        return(NULL)
    }
    
    n <- x
    if(!isScalarNumeric(x)){
        if(is.numeric(x)){
            set <- tryCatch({
                    colorRampPalette(set)
                    set
                }, error = function(e){
                    warning("Please specify correct color for ", 
                            "colorRampPalette. Using default. ",
                            "\n  Error was:\n    ", e)
                    NULL
            })
            return(set)
        }
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
    
    if(length(x) > 0){
        ans <- cols[x]
        names(ans) <- as.character(x)
        return(ans)
    }
    
    return(cols)
}

#' 
#' Creates a color mapping for the given annotation or return NULL
#' 
#' @noRd
getAnnoColors <- function(colorSet, annotation){
    if(!is.character(colorSet) & is.null(annotation)){
        return(NULL)
    }
    mapply(grp=annotation, set=colorSet, SIMPLIFY=FALSE, 
            FUN=function(grp, set) getXColors(grp[!duplicated(grp)], set))
}

#'
#' Ignore warning of the text aestethics in ggplot for plotly since it is 
#' inofficial till now
#' 
#' @noRd
ignoreAesTextWarning <- function(expr){
    withCallingHandlers(expr, 
            warning = function(w){
                if(endsWith(conditionMessage(w), "unknown aesthetics: text")){
                    invokeRestart("muffleWarning")
                }
    })
}

#' 
#' Get annotation for row/col dependent on input
#' @noRd
getAnnoHeatmap <- function(x, matrix, groups, nClust, extractFun=colData){
    
    # select annotations based on metadata (colData/mcols)
    if(!isScalarNA(groups)){
        ans <- as.data.frame(extractFun(x)[, groups])
        colnames(ans) <- groups
    } else {
        ans <- as.data.frame(extractFun(x)[,character()])
    }
    
    # add clustering
    if(isScalarNumeric(nClust) & isScalarNA(groups)){
        clusters <- cutree(hclust(dist(matrix)), nClust)
        ans[["nClust"]] <- as.character(clusters)
    }
    
    # return NULL if no annotation is requested
    if(ncol(ans) == 0){
        return(NULL)
    }
    
    # due to not propagating rownames correctly in mcols/rowData in R < 3.5.0
    # check if rownames is null and then retrive from object
    rownames(ans) <- rownames(extractFun(x))
    if(is.null(rownames(extractFun(x)))){
        if(ncol(x) == nrow(ans)){
            rownames(ans) <- colnames(x)
        } else {
            rownames(ans) <- rownames(x)
        }
    }
    return(ans)
}

getNiceName <- function(x, maxChar=12){
    stopifnot(maxChar > 2)
    ifelse(nchar(x) > maxChar, 
            paste0(substr(x, 0, maxChar - 2), ".."),
            x)
}

#' 
#' This function is used by the plotDispEsts function.
#' 
#' TODO
#' 
#' @noRd
getDispEstsData <- function(ods, mu=NULL){
    if(is.null(theta(ods))){
        stop('Please fit the ods first. ods <- fit(ods)')
    }
    odsMu <- rowMeans(counts(ods, normalized=TRUE))
    if(is.null(mu)){
        mu <- odsMu
    }
    theta <- theta(ods)
    xidx <- 10^(seq.int(max(-5,log10(min(mu))-1), log10(max(mu))+0.1, 
            length.out = 500))
    
    # fit DESeq2 parametric Disp Fit
    fit <- parametricDispersionFit(mu, 1/theta)
    pred <- fit(xidx)
    return(list(
        mu=mu,
        disp=theta,
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

#' @rdname getter_setter_functions
#' @export
getBestQ <- function(ods){
    if('optimalEncDim' %in% names(metadata(ods))){
        return(metadata(ods)[['optimalEncDim']])
    }
    
    if('encDimTable' %in% names(metadata(ods))){
        encTable <- metadata(ods)[['encDimTable']]
        return(getBestQDT(encTable, 'aucPR'))
    }
    # warning('Please find the optimal encoding dimension by running. ')
    return(NA_integer_)
}

getBestQDT <- function(dt, usedEvalMethod='aucPR', digits=10){
    if('evalMethod' %in% colnames(dt)){
        testFun <- ifelse(all(dt[,evalMethod == usedEvalMethod]), 
                which.max, which.min)
    } else {
        testFun <- which.max
    }
    
    dt[,encodingDimension[
            seq_len(.N) == testFun(round(evaluationLoss, digits))]]
}

#'
#' Estimation of Q
#' 
#' Estimating the best q for the given data set
#' 
#' @param ods An OutriderDataSet object
#' @return The estimated dimension of hidden confounders
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' 
#' estimateBestQ(ods)
#' 
#' @export
estimateBestQ <- function(ods){
    round(max(2, min(500, 3.7 + 0.16*ncol(ods))))
}
