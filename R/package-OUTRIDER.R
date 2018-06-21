#'
#' This is the import mapping for the scared package
#' 
#' @noRd
#' 
#' @name OUTRIDER
#'
#' @import data.table
#' 
#' @importFrom Biobase rowMax
#' 
#' @importFrom BiocGenerics estimateSizeFactors plotDispEsts
#' 
#' @importFrom DESeq2 normalizationFactors normalizationFactors<-
#'          sizeFactors sizeFactors<- counts counts<-
#'          DESeqDataSetFromMatrix DESeqDataSet
#'          makeExampleDESeqDataSet show fpkm fpm
#' 
#' @importFrom SummarizedExperiment colData colData<- assays assays<- 
#'          assayNames mcols mcols<- width rowData assay assay<-
#' 
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' 
#' @importFrom GenomicRanges reduce
#' 
#' @importFrom BiocParallel bplapply bpparam MulticoreParam
#' 
#' @importFrom tidyr %>%
#' 
#' @importFrom stats p.adjust setNames sd dnbinom quantile optim var pnbinom 
#'          dnbinom median ppoints qbeta runif cor cutree hclust dist lm predict
#'          rnorm glm Gamma
#' 
#' @importFrom methods validObject is new as 
#' 
#' @importFrom ggplot2 ggplot geom_point aes geom_hline geom_vline 
#'          geom_vline ggtitle theme element_text geom_tile
#'          geom_segment scale_colour_manual scale_y_log10 labs unit
#'          geom_histogram scale_x_log10 scale_fill_manual
#' 
#' @importFrom scales trans_format math_format         
#'         
#' @importFrom DelayedArray rowMins
#' 
#' @importFrom S4Vectors DataFrame rowSums colSums metadata metadata<-
#' 
#' @importFrom utils read.table head compareVersion
#' 
#' @importFrom graphics plot abline grid points mtext polygon axis text box
#'         legend title lines par
#' 
#' @importFrom gplots barplot2 heatmap.2
#' 
#' @importFrom matrixStats rowSds rowMedians rowQuantiles
#' 
#' @importFrom plotly plot_ly add_markers ggplotly plotly_build style layout
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @importFrom BBmisc isScalarLogical isScalarNumeric isScalarInteger 
#'         isScalarCharacter isScalarValue isFALSE
#' 
#' @importFrom reticulate import py_get_attr
#' 
#' @importFrom LSD heatscatter
#' 
#' @importFrom pcaMethods pca loadings
#' 
#' @importFrom compiler cmpfun
#' 
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib OUTRIDER
#' 
# ###
# # TODO old package to be checked again!
#' 
#' @importFrom gplots  bluered
#' @importFrom tidyr gather
#' @importFrom stats as.formula complete.cases coefficients relevel
#' @importFrom plyr join .
#' 
#' 

NULL

#' 
#' TODO This is to get rid of the warnings of undefined variables
#'     due to the nature of data.table and ggplot/plotly
#' @noRd
globalVariables(c("disp", "padjust", "pValue", "zScore", "padj_rank", "mu", 
        "l2fc", "sampleID", "normcounts", "sample.rank", "N2", "color", "obs",
        "Var1", "Var2", "value", "medianCts", "norm_rank", ".x",
        "passedFilter", "loggeomeans"), package="OUTRIDER")

