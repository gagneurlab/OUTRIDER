#'
#' This is the import mapping for the scared package
#' 
#' @noRd
#' 
#' @name OUTRIDER
#'
#' @import data.table
#' 
#' @importFrom methods as is new validObject
#' 
#' @importFrom Biobase rowMax
#' 
#' @importFrom BiocGenerics estimateSizeFactors plotDispEsts
#' 
#' @importFrom DESeq2 normalizationFactors normalizationFactors<-
#'          sizeFactors sizeFactors<- counts counts<-
#'          DESeqDataSetFromMatrix DESeqDataSet
#'          makeExampleDESeqDataSet show fpkm fpm
#'          estimateSizeFactorsForMatrix replaceOutliers
#'          dispersions
#' 
#' @importFrom SummarizedExperiment colData colData<- assays assays<-
#'          assayNames mcols mcols<- assay assay<-
#' 
#' @importFrom BBmisc isScalarLogical isScalarNumeric isScalarCharacter isFALSE
#'          isScalarValue chunk
#' 
#' @importFrom BiocParallel bplapply bpparam SerialParam bpisup bpstart bpstop
#' 
#' @importFrom compiler cmpfun
#' 
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' 
#' @importFrom GenomicRanges GRanges reduce width 
#' 
#' @importFrom ggplot2 aes geom_histogram geom_smooth geom_tile geom_point 
#'          ggplot labs scale_x_log10 scale_fill_manual theme ylim ggtitle
#'          geom_vline geom_text
#' 
#' @importFrom ggpubr grids
#'          
#' @importFrom glmnet cv.glmnet
#' 
#' @importFrom lbfgs lbfgs
#'          
#' @importFrom gplots barplot2 bluered heatmap.2
#' 
#' @importFrom graphics plot abline axis box grid legend lines mtext par points 
#'          polygon text title
#'
#' @importFrom IRanges IRanges
#' 
#' @importFrom matrixStats rowSds rowMedians rowQuantiles rowMeans2
#' 
#' @importFrom pcaMethods pca loadings
#' 
#' @importFrom plotly %>% ggplotly layout plotly_build plot_ly
#' 
#' @importFrom plyr .
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @importFrom Rcpp sourceCpp
#' 
#' @importFrom reticulate import py_get_attr
#' 
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' 
#' @importFrom scales math_format trans_format
#' 
#' @importFrom splines bs
#' 
#' @importFrom stats cor coefficients cutree dist dnbinom hclust p.adjust 
#'          setNames sd optimize rlnorm
#'          quantile optim var pnbinom  median ppoints qbeta runif 
#'            lm predict rnorm glm Gamma rnbinom
#' 
#' @importFrom utils read.table head compareVersion
#' 
#' @importFrom PRROC pr.curve
#' 
#' @useDynLib OUTRIDER
#' 
NULL

#' 
#' TODO This is to get rid of the warnings of undefined variables
#'     due to the nature of data.table and ggplot/plotly
#' @noRd
NULL
globalVariables(c(
    "enc",
    "eva",
    "encodingDimension", 
    "evalMethod", 
    "evaluationLoss",
    "normcounts",
    "opt",
    "padj_rank", 
    "padjust", 
                  
## color", "frac", "Fraction", "ExprType", 
##       "evaluationLoss", "loggeomeans", "lty", 
##        "lwd", "medianCts", "mu", "negLog10pVal", "normcounts", "norm_rank", 
##        "obs", "onlyFull","passedFilter", "pValue", 
## "zScore"), 
       "sampleID", "value", "V1", "Var1", "Var2", ".x", 
        "z"), package="OUTRIDER")

