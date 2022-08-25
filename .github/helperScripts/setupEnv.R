BTYPE <- ifelse(.Platform$OS.type == 'unix', "source", "both")
NCPUS <- 3
Sys.setenv(MAKEFLAGS = "-j3")
BIOC_VERSION <- Sys.getenv("BIOC_VERSION")
START_TIME <- Sys.time()

print_log <- function(...){
    hash_line <- paste0(rep("#", 10), collapse="")
    message(paste0("\n", hash_line, "\n### ", date(), ": ", ..., "\n", hash_line, "\n"))
}

# install Bioconductor
if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", Ncpus=NCPUS)
BiocManager::install("BiocVersion", version=BIOC_VERSION,ask = FALSE)

INSTALL <- BiocManager::install
installIfReq <- function(p, Ncpus=NCPUS, ...){
    for(j in p)
        if(!requireNamespace(j, quietly=TRUE))
            BiocManager::install(j, Ncpus=Ncpus, ...)
}

# since the current XML package is not compatible with 3.6 anymore
if(!requireNamespace("XML", quietly=TRUE) & R.version[['major']] == "3"){
    installIfReq(p="devtools", type=BTYPE, Ncpus=NCPUS)
    devtools::install_version("XML", version="3.99-0.3")
}

# install needed packages
# add testthat to pre installation dependencies due to: https://github.com/r-lib/pkgload/issues/89
for(p in c("textshaping","Biostrings","XML", "xml2", "testthat", "devtools", "covr", "roxygen2", "BiocCheck", "R.utils", 
            "GenomeInfoDbData", "rtracklayer", "hms")){
    installIfReq(p=p, type=BTYPE, Ncpus=NCPUS)
}

# because of https://github.com/r-windows/rtools-installer/issues/3
if("windows" == .Platform$OS.type){
    print_log("Install XML on windows ...")
    BTYPE <- "win.binary"
    installIfReq(p=c("XML", "xml2", "RSQLite", "progress", "AnnotationDbi", "BiocCheck"))
    
    print_log("Install source packages only for windows ...")
    INSTALL(c("GenomeInfoDbData", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene"),
            type="both", version=BIOC_VERSION)                                                                          

} else {
    BTYPE <- "source"
}

print_log("Update Packages")
BiocManager::install(ask=FALSE, Ncpus=NCPUS, version=BIOC_VERSION)

print_log("Install dev package")
devtools::install(".", dependencies=TRUE)


# fix knitr for 3.6 for more details see BiocStyle issue 78
# https://github.com/Bioconductor/BiocStyle/issues/78
if(R.version[['major']] == "3"){
    options(repos=c(CRAN="http://cran.rstudio.com"))
    installIfReq(p="devtools", type=BTYPE, Ncpus=NCPUS)
    devtools::install_version("knitr", version="1.29", type="source")
}

# to get OUTRIDER session info
try({ library(OUTRIDER) })
print(BiocManager::valid())
