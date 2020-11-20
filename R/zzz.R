
# global reference to py_outrider (will be initialized in .onLoad)
py_outrider <- NULL

#'
#' The .onLoad function which is run during package loading
#' 
#' @noRd
#' 
.onLoad <- function(libname, pkgname){
    path <- system.file("python", package = packageName())
    py_outrider <<- reticulate::import_from_path(module="py_outrider",
                                                path=path) #, delay_load=TRUE)
    # py_outrider <<- reticulate::import("py_outrider", delay_load=TRUE)
}