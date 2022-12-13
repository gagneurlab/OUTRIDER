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

installIfReq <- function(p, Ncpus=NCPUS, ...){
    for(j in p)
        if(!requireNamespace(j, quietly=TRUE))
            BiocManager::install(j, Ncpus=Ncpus, ...)
}

# install needed packages
# add testthat to pre installation dependencies due to: https://github.com/r-lib/pkgload/issues/89
for(p in c("Biostrings","XML", "xml2", "testthat", "devtools", "covr", "roxygen2", "BiocCheck", "R.utils", 
            "GenomeInfoDbData", "rtracklayer", "hms", "RSQLite", "progress", "AnnotationDbi")){
    installIfReq(p=p, Ncpus=NCPUS)
}

# because of https://github.com/r-windows/rtools-installer/issues/3
if("windows" == .Platform$OS.type){
    print_log("Install source packages only for windows ...")
    BiocManager::install(c("GenomeInfoDbData", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene"),
            type="both", version=BIOC_VERSION)                                                                          
}

print_log("Update Packages")
BiocManager::install(ask=FALSE, Ncpus=NCPUS, version=BIOC_VERSION)

print_log("Install dev package")
devtools::install(".", dependencies=TRUE)

# to get OUTRIDER session info
try({ library(OUTRIDER) })
print(BiocManager::valid())
