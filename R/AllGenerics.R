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



### OUTRIDER2 generics

#' @rdname getter_setter_functions
#' @export
setGeneric("modelParams", function(object, paramName){
    standardGeneric("modelParams")
}) 


#' @rdname getter_setter_functions
#' @export
setGeneric("modelParams<-", signature = "object", 
           function(object, value, paramName) standardGeneric("modelParams<-"))


#' @rdname getter_setter_functions
#' @export
setGeneric("observed", function(object, normalized=FALSE, minE=0, ...){ 
    standardGeneric("observed")
})

#' @rdname getter_setter_functions
#' @export
setGeneric("observed<-", signature = "object", 
           function(object, value, ...) standardGeneric("observed<-"))


#' @rdname getter_setter_functions
#' @export
setGeneric("expected", function(object, ...){ 
    standardGeneric("expected")
})

#' @rdname getter_setter_functions
#' @export
setGeneric("expected<-", signature = "object", 
           function(object, value, ...) standardGeneric("expected<-"))


#' @rdname getter_setter_functions
#' @export
setGeneric("preprocessed", function(object, ...) ("preprocessed"))

#' @rdname getter_setter_functions
#' @export
setGeneric("preprocessed<-", signature = "object", 
           function(object, value, ...) standardGeneric("preprocessed<-"))
