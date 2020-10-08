BTYPE <- ifelse(.Platform$OS.type == 'unix', "source", "both")
NCPUS <- 6
START_TIME <- Sys.time()

print_log <- function(...){
    hash_line <- paste0(rep("#", 10), collapse="")
    message(paste0("\n", hash_line, "\n### ", date(), ": ", ..., "\n", hash_line, "\n"))
}

installIfReq <- function(p, type=BTYPE, Ncpus=NCPUS, ...){
    for(j in p){
        if(!requireNamespace(j, quietly=TRUE)){
            print_log("Install ", j)
            INSTALL(j, type=type, Ncpus=Ncpus, ...)
        }
    }
}

# install Bioconductor
if(!requireNamespace("BiocManager", quietly=TRUE)){
    print_log("Install BiocManager")
    install.packages("BiocManager", Ncpus=NCPUS)
}
INSTALL <- BiocManager::install

# since the current XML package is not compatible with 3.6 anymore
if(!requireNamespace("XML", quietly=TRUE) & R.version[['major']] == "3"){
    installIfReq(p="devtools", type=BTYPE, Ncpus=NCPUS)
    devtools::install_version("XML", version="3.99-0.3")
}

# because of https://github.com/r-windows/rtools-installer/issues/3
if("windows" == .Platform$OS.type){
    print_log("Install XML on windows ...")
    BTYPE <- "win.binary"
    installIfReq(p=c("XML", "xml2", "RSQLite", "progress", "AnnotationDbi", "BiocCheck"))
    
    print_log("Install source packages only for windows ...")
    INSTALL(c("GenomeInfoDbData", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene"), type="both")
} else {
    BTYPE <- "source"
}

# install needed packages
# add testthat to pre installation dependencies due to: https://github.com/r-lib/pkgload/issues/89
for(p in c("XML", "xml2", "testthat", "devtools", "covr", "roxygen2", "BiocCheck", "R.utils", 
            "GenomeInfoDbData", "rtracklayer", "hms")){
    installIfReq(p=p, type=BTYPE, Ncpus=NCPUS)
}

# install OUTRIDER with its dependencies with a timeout due to 
# travis (50 min) and appveyor (60 min) set installation warmup to 30 min max
maxTime <- max(30, (60*30 - difftime(Sys.time(), START_TIME, units="sec")))
R.utils::withTimeout(timeout=maxTime, {
    try({
        print_log("Update packages")
        INSTALL(ask=FALSE, type=BTYPE, Ncpus=NCPUS)
    
        print_log("Install OUTRIDER")
        devtools::install(".", dependencies=TRUE, upgrade=TRUE, 
                type=BTYPE, Ncpus=NCPUS)
    })
})

# fix knitr for 3.6 for more details see BiocStyle issue 78
# https://github.com/Bioconductor/BiocStyle/issues/78
if(R.version[['major']] == "3"){
    installIfReq(p="devtools", type=BTYPE, Ncpus=NCPUS)
    devtools::install_version("knitr", version="1.29", type="source")
}
