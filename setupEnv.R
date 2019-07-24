# install Bioconductor dependent on the R version
if(0 < compareVersion("3.5.0", paste(R.Version()[c("major", "minor")], collapse="."))){
    if(isFALSE(requireNamespace("BiocInstaller"))){
        source("https://bioconductor.org/biocLite.R")
        biocLite(c("BiocInstaller"), Ncpus=6)
    }
    INSTALL <- BiocInstaller::biocLite
} else {
    if(isFALSE(requireNamespace("BiocInstaller")))
        install.packages("BiocManager", Ncpus=6)
    INSTALL <- BiocManager::install
}

# install needed packages
for(p in c("devtools", "covr", "roxygen2", "BiocCheck", "R.utils"))
    if(isFALSE(requireNamespace(p, quietly=TRUE)))
        INSTALL(p, Ncpus=6)

# install OUTRIDER with its dependencies with a timeout due to travis 50 min
R.utils::withTimeout(timeout=2400, {
    INSTALL(ask=FALSE, Ncpus=6)
    devtools::install(".", ask=FALSE, dependencies=TRUE, upgrade=TRUE)
})

