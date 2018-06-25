
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


#' 
#' Read/Write NB Model
#' 
#' Writes out to a file or reads in from a file a fitted NB model.
#' The model contains:
#' \itemize{
#'     \item{feature name}
#'     \item{log geometric means}
#'     \item{fitted mean}
#'     \item{dispersion}
#' }
#' 
#' 
#' @param ods OutriderDataSet
#' @param modelFile specify file, where model should be saved.
#' @return An OutriderDataSet containing the loaded model or None if its saved.
#' 
#' @examples 
#' # fit and write out model
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- fit(ods)
#' writeNBModel(ods, modelFile="model.tsv")
#' 
#' # read saved model
#' ods <- makeExampleOutriderDataSet()
#' ods <- readNBModel(ods, 'model.tsv')
#' ods <- estimateSizeFactors(ods)
#' ods <- computePvalues(ods)
#' 
#' @rdname readNBModel
#' @aliases writeNBModel
#' @export
writeNBModel <- function(ods, modelFile){
    if(!isScalarCharacter(modelFile, na.ok=FALSE)){
        stop("Please provide a correct file name (Character) for the model.")
    }
    if(is.null(rownames(ods))){
        stop("Please provide rownames, if you want to save the fit model.")
    }
    if(any(duplicated(rownames(ods)))){
        stop("Please provide uniq rownames, if you want to save the fit model.")
    }
    if(!dir.exists(dirname(modelFile))){
        dir.create(dirname(modelFile), recursive=TRUE)
    }
    cols2save <- c('loggeomeans', 'mu', 'disp')
    if(any(!cols2save %in% colnames(mcols(ods)))){
        missedCols <- paste(
                cols2save[!cols2save %in% colnames(mcols(ods))], collapse=", ")
        stop(paste0("The following mcols are missing: ", missedCols,
                ". Please run first estimateSizeFactors() and fit()."))
    }
    
    model2save <- cbind(
        data.table(rownames=rownames(ods)),
        as.data.table(mcols(ods)[cols2save]))
    fwrite(model2save, modelFile, sep="\t")
    return(invisible())
}

#' @rdname readNBModel
#' @export
readNBModel <- function(ods, modelFile){
    if(!isScalarCharacter(modelFile)){
        stop("Please provide a correct file name (Character) for the model.")
    }
    if(!file.exists(modelFile)){
        stop(paste0("Can not find given model file: '", modelFile, "'."))
    }
    if(is.null(rownames(ods))){
        stop("Please provide rownames to be able to match the model")
    }
    model <- fread(modelFile)
    if(!all(rownames(ods) %in% model[,rownames])){
        failedFeatures <- rownames(ods)[rownames(ods) %in% model[,rownames]]
        stop("We could not find all needed fit parameters in the model. ",
                "We are missing the following features: ", 
                paste(collapse="\n", head(failedFeatures, 20)))
    }
    
    # overwrite fit parameters of object
    setkey(model, rownames)
    mcols(ods)[c("loggeomeans", "mu", "disp")] <- DataFrame(
            model[rownames(ods), .(loggeomeans, mu, disp)])
    validObject(ods)
    
    return(ods)
}


#'
#' This function is used by the plotDispEsts function.
#' 
#' TODO
#' 
#' @noRd
getDispEstsData <- function(ods, mu=NULL){
    odsMu <- rowMeans(counts(ods, normalized=TRUE))
    if(is.null(mu)){
        mu <- odsMu
    }
    disp <- mcols(ods)$disp
    xidx <- 10^(seq.int(log10(min(mu))-1, log10(max(mu))+0.1, length.out = 500))
    
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
