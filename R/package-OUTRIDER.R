#'
#' This is the import mapping for the OUTRIDER package
#' 
#' @noRd
#' 
#' @name OUTRIDER
#'
#' @rawNamespace import(data.table, except=melt)
#' 
#' @import methods
#' 
#' @importFrom BiocGenerics estimateSizeFactors plotDispEsts
#' 
#' @importFrom DESeq2 normalizationFactors normalizationFactors<-
#'          sizeFactors sizeFactors<- counts counts<-
#'          DESeqDataSetFromMatrix DESeqDataSet
#'          makeExampleDESeqDataSet show fpkm fpm
#'          estimateSizeFactorsForMatrix replaceOutliers dispersions
#' 
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#'          assays assays<- assayNames mcols mcols<- assay assay<- rowRanges
#' 
#' @importFrom BBmisc isScalarLogical isScalarNumeric isScalarCharacter isFALSE
#'          isScalarValue isScalarNA chunk seq_col seq_row
#' 
#' @importFrom BiocParallel bplapply bpparam SerialParam bpisup bpstart bpstop
#'          bpmapply
#' 
#' @importFrom GenomicFeatures exonsBy
#'
#' @importFrom txdbmaker makeTxDbFromGFF
#' 
#' @importFrom GenomicRanges GRanges reduce width start end findOverlaps 
#' 
#' @importFrom ggplot2 ggplot aes annotate geom_bar geom_histogram 
#'          geom_hline geom_smooth geom_point labs scale_x_log10 
#'          scale_y_log10 scale_fill_manual scale_color_manual 
#'          scale_fill_brewer scale_color_brewer theme ylim 
#'          ggtitle geom_vline geom_text scale_linetype_manual geom_line 
#'          geom_abline theme_bw element_blank xlab ylab scale_color_identity 
#'          facet_grid facet_wrap quo_name
#' 
#' @importFrom ggrepel geom_text_repel
#' 
#' @importFrom grDevices colorRampPalette
#' 
#' @importFrom pheatmap pheatmap
#' 
#' @importFrom heatmaply heatmaply
#' 
#' @importFrom graphics plot abline axis box grid legend lines mtext par points 
#'          polygon text title
#'
#' @importFrom IRanges IRanges
#' 
#' @importFrom matrixStats rowSds rowMedians rowQuantiles rowMeans2 rowCumsums
#' 
#' @importFrom pcaMethods pca loadings
#' 
#' @importFrom plotly %>% ggplotly layout plotly_build plot_ly
#' 
#' @importFrom plyr .
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @importFrom reshape2 melt
#' 
#' @importFrom S4Vectors DataFrame metadata metadata<- endoapply queryHits
#'          subjectHits values values<-
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
#' @importFrom utils read.table head compareVersion packageVersion
#' 
#' @importFrom PRROC pr.curve
#' 
#' @useDynLib OUTRIDER
#' 
NULL

#' 
#' TODO This is to get rid of the warnings of undefined variables
#'     due to the nature of data.table and ggplot/plotly
#' For reference:
#' https://www.r-bloggers.com/no-visible-binding-for-global-variable/
#' @noRd
globalVariables(package="OUTRIDER", c(
        "color",
        "enc",
        "encodingDimension", 
        "eva",
        "evalMethod", 
        "evaluationLoss",
        "expected",
        "expRank",
        "ExprType", 
        "feature_id",
        "fill",
        "frac",
        "Fraction",
        "GENE_ID",
        "group",
        "lty", 
        "lwd", 
        "medianCts", 
        "negLog10pVal",
        "normcounts",
        "normCts",
        "norm_rank", 
        "opt",
        "padj_rank", 
        "padjust", 
        "obs",
        "observed",
        "onlyFull",
        "passedFilter",
        "Rank",
        "rawcounts",
        "sampleID",
        "sizeFactor",
        "value",
        "variable",
        "V1",
        "Var1",
        "Var2",
        ".x",
        "y",
        "z",
        "AberrantByGene",
        "AberrantBySample",
        "FDR_set",
        "aberrantLabel",
        "geneID"))

