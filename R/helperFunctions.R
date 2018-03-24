
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
        cols <- cols[1:n]
    }
    
    if(length(x) > 1){
        ans <- cols[x]
        names(ans) <- as.character(x)
        return(ans)
    }
    
    return(cols)
}

#' 
#' heatmap.2 without a trace as default.
#' 
#' @examples 
#'   heatmapNotrace(matrix(runif(100), nrow=10))
#'   
#' @author baderd
#' @noRd
heatmapNotrace <- function( x, denscol='green', col=bluered(50), 
                            key.par=list(las=1), keysize=1, ...){
    requireNamespace("gplots")
    heatmap.2(
        x=x, trace='none', 
        key.ylab='', key.title='', 
        keysize=keysize, key.par=key.par,
        col=col, denscol=denscol
        ,...
    )
}


#' 
#' Read/Write NB Model
#' 
#' Writes out to a file or reads in from a file a fitted NB model.
#' The model contains:
#'   * feature name
#'   * log geometric means
#'   * fitted mean
#'   * dispersion
#' 
#' @param ods OutriderDataSet
#' @param modelFile specify file, where model should be saved.
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
    
    model2save <- cbind(
        data.table(rownames=rownames(ods)),
        as.data.table(rowData(ods)[c('loggeomeans', 'mu', 'disp')]))
    fwrite(model2save, modelFile, sep="\t")
    return()
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
