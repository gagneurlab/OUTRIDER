#' @rdname aberrant
#' @export
setGeneric("aberrant", function(object, ...) standardGeneric("aberrant"))

#' @importFrom DESeq2 dispersions
#' @noRd
#' @export
DESeq2::dispersions

#' @rdname filterExpression
#' @export
setGeneric("filterExpression", function(object, ...)
        standardGeneric("filterExpression"))

#' @importFrom generics fit
#' @noRd
#' @export
generics::fit

#' @rdname plotFunctions
#' @export
setGeneric("plotAberrantPerSample", function(object, ...)
        standardGeneric("plotAberrantPerSample"))

#' @rdname plotFunctions
#' @export
setGeneric("plotCountCorHeatmap", function(object, ...) 
        standardGeneric("plotCountCorHeatmap"))

#' @rdname plotFunctions
#' @export
setGeneric("plotManhattan", function(object, ...) 
    standardGeneric("plotManhattan"))

#' @importFrom DESeq2 plotDispEsts
#' @noRd
#' @export
DESeq2::plotDispEsts

#' @rdname plotFunctions
#' @export
setGeneric("plotEncDimSearch", function(object, ...)
        standardGeneric("plotEncDimSearch"))

#' @rdname plotFunctions
#' @export
setGeneric("plotQQ", function(object, ...) standardGeneric("plotQQ"))

#' @rdname plotFunctions
#' @export
setGeneric("plotVolcano", function(object, ...) standardGeneric("plotVolcano"))

#' @rdname results
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))
