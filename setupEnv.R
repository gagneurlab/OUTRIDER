print_log <- function(...){
    hash_line <- paste0(rep("#", 10), collapse="")
    message(paste0("\n", hash_line, "\n### ", date(), ": ", ..., "\n", hash_line, "\n"))
}

# install Bioconductor dependent on the R version
R_VERSION <- paste(R.Version()[c("major", "minor")], collapse=".")
print_log("Current R version: ", R_VERSION)
if(0 < compareVersion("3.5.0", R_VERSION)){
    if(!requireNamespace("BiocInstaller", quietly=TRUE)){
        print_log("Install BiocInstaller")
        source("https://bioconductor.org/biocLite.R")
        biocLite(c("BiocInstaller"), Ncpus=6)
    }
    INSTALL <- BiocInstaller::biocLite
} else {
    if(!requireNamespace("BiocManager", quietly=TRUE)){
        print_log("Install BiocManager")
        install.packages("BiocManager", Ncpus=6)
    }
    
    INSTALL <- BiocManager::install
}

# install needed packages
# add testthat to pre installation dependencies due to: https://github.com/r-lib/pkgload/issues/89
for(p in c("testthat", "devtools", "covr", "roxygen2", "BiocCheck", "R.utils")){
    if(!requireNamespace(p, quietly=TRUE)){
        print_log("Install ", p)
        INSTALL(p, Ncpus=6)
    }
}

# install OUTRIDER with its dependencies with a timeout due to travis 50 min
R.utils::withTimeout(timeout=2400, {
    print_log("Update packages")
    INSTALL(ask=FALSE, Ncpus=6)
    
    print_log("Install OUTRIDER")
    devtools::install(".", ask=FALSE, dependencies=TRUE, upgrade=TRUE, Ncpus=6)
})
