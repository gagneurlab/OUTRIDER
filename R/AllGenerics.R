#'
#' OUTRIDER specific functions defined as generics
#' 
#' @noRd
NULL


#' @rdname results
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))


#' @rdname getter_setter_functions
#' @export
setGeneric("modelParams", function(object, paramName) standardGeneric("modelParams"))

#' @rdname getter_setter_functions
#' @export
setGeneric("modelParams<-", signature = "object", 
           function(object, paramName, value) standardGeneric("modelParams<-"))


#' @rdname getter_setter_functions
#' @export
setGeneric("observed", function(object, normalized) standardGeneric("observed"))

#' @rdname getter_setter_functions
#' @export
setGeneric("observed<-", signature = "object", 
           function(object, value) standardGeneric("observed<-"))

#' @rdname getter_setter_functions
#' @export
setGeneric("preprocessed", function(object) ("preprocessed"))

#' @rdname getter_setter_functions
#' @export
setGeneric("preprocessed<-", signature = "object", 
           function(object, value) standardGeneric("preprocessed<-"))