
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
    
    cols <- suppressWarnings(brewer.pal(n=n, name=set))
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
        stop(paste("We could not find all needed fit parameters in the model.",
                "We are missing the following features: ", paste(collapse="\n",
                    head(rownames(ods)[rownames(ods) %in% model[,rownames]], 20)
                )))
    }
    
    # overwrite fit parameters of object
    setkey(model, rownames)
    mcols(ods)[c("loggeomeans", "mu", "disp")] <- DataFrame(
            model[rownames(ods), .(loggeomeans, mu, disp)])
    validObject(ods)
    
    return(ods)
}
