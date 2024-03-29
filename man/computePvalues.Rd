% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pValMatrix.R
\docType{methods}
\name{computePvalues}
\alias{computePvalues}
\alias{computePvalues,OutriderDataSet-method}
\title{Calculate P-values}
\usage{
computePvalues(object, ...)

\S4method{computePvalues}{OutriderDataSet}(
  object,
  alternative = c("two.sided", "greater", "less"),
  method = "BY",
  subsets = NULL,
  BPPARAM = bpparam()
)
}
\arguments{
\item{object}{An OutriderDataSet}

\item{...}{additional params, currently not used.}

\item{alternative}{Can be one of "two.sided", "greater" or "less" to specify
the alternative hypothesis used to calculate the P-values,
defaults to "two.sided"}

\item{method}{Method used for multiple testing}

\item{subsets}{A named list of named lists specifying any number of gene 
subsets (can differ per sample). For each subset, FDR correction
will be limited to genes in the subset, and the FDR corrected 
pvalues stored as an assay in the ods object in addition to the 
transcriptome-wide FDR corrected pvalues. See the examples for 
how to use this argument.}

\item{BPPARAM}{Can be used to parallelize the computation, defaults to
bpparam()}
}
\value{
An OutriderDataSet object with computed nominal and adjusted P-values
}
\description{
This function computes the P-values based on the fitted negative binomial
model being an outlier for the given sample gene combination. It computes
two matrices with the same dimension as the count matrix (samples x genes),
which contain the corresponding P-values and adjusted P-values of every 
count. The adjusted P-values are computed across all genes per sample.
}
\examples{
ods <- makeExampleOutriderDataSet()
\dontshow{
    ods <- ods[1:10,1:10]
}
ods <- estimateSizeFactors(ods)
ods <- fit(ods)

ods <- computePvalues(ods)

assays(ods)[['pValue']][1:10,1:10]

# example of restricting FDR correction to subsets of genes of interest
genesOfInterest <- list("sample_1"=sample(rownames(ods), 3), 
                         "sample_2"=sample(rownames(ods), 8), 
                         "sample_6"=sample(rownames(ods), 5))
ods <- computePvalues(ods, subsets=list("exampleSubset"=genesOfInterest))
padj(ods, subsetName="exampleSubset")[1:10,1:10]
ods <- computePvalues(ods, 
                 subsets=list("anotherExampleSubset"=rownames(ods)[1:5]))
padj(ods, subsetName="anotherExampleSubset")[1:10,1:10]

}
\seealso{
p.adjust
}
