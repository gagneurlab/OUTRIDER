# install Bioconductor dependent on the R version
R_VERSION <- paste(R.Version()[c("major", "minor")], collapse="."))
message(paste0(date(), ": Current R version: ", R_VERSION))
if(0 < compareVersion("3.5.0", R_VERSION)){
    if(!requireNamespace("BiocInstaller")){
        message("Install BiocInstaller")
        source("https://bioconductor.org/biocLite.R")
        biocLite(c("BiocInstaller"), Ncpus=6)
    }
    INSTALL <- BiocInstaller::biocLite
} else {
    if(!requireNamespace("BiocInstaller")){
        message("Install BiocManager")
        install.packages("BiocManager", Ncpus=6)
    }
    INSTALL <- BiocManager::install
}

# install needed packages
for(p in c("devtools", "covr", "roxygen2", "BiocCheck", "R.utils")){
    if(!requireNamespace(p, quietly=TRUE)){
        message(paste("Install", p))
        INSTALL(p, Ncpus=6)
    }
}

# install OUTRIDER with its dependencies with a timeout due to travis 50 min
R.utils::withTimeout(timeout=2400, {
    message("Update packages")
    INSTALL(ask=FALSE, Ncpus=6)
    
    message("Install OUTRIDER")
    devtools::install(".", ask=FALSE, dependencies=TRUE, upgrade=TRUE)
})

