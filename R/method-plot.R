#'
#' @title Visualization functions for OUTRIDER
#'
#' @description The OUTRIDER package provides mutliple functions to visualize
#' the data and the results of a full data set analysis.
#'
#' This is the list of all plotting function provided by OUTRIDER:
#' \itemize{
#'   \item plotAberrantPerSample()
#'   \item plotVolcano()
#'   \item plotExpressionRank()
#'   \item plotQQ()
#'   \item plotExpectedVsObservedCounts()
#'   \item plotCountCorHeatmap()
#'   \item plotCountGeneSampleHeatmap()
#'   \item plotSizeFactors()
#'   \item plotFPKM()
#'   \item plotExpressedGenes()
#'   \item plotDispEsts()
#'   \item plotPowerAnalysis()
#'   \item plotEncDimSearch()
#' }
#'
#' For a detailed description of each plot function please see the details.
#' Most of the functions share the same parameters.
#'
#### Data specific parameters
#' @param ods,object An OutriderDataSet object.
#' @param sampleID,geneID A sample or gene ID, which should be plotted.
#'             Can also be a vector. Integers are treated as indices.
#' @param padjCutoff,zScoreCutoff Significance or Z-score cutoff
#'             to mark outliers
#' @param global Flag to plot a global Q-Q plot, default FALSE
#' @param outlierRatio The fraction to be used for the outlier sample filtering
#' @param normalized If TRUE, the normalized counts are used, the default,
#'             otherwise the raw counts
#' @param compareDisp If TRUE, the default, and if the autoCorrect normalization
#'             was used it computes the dispersion without autoCorrect and
#'             plots it for comparison.
#### Graphical parameters
#' @param main Title for the plot, if missing a default title will be used.
#' @param groups A character vector containing either group assignments of
#'             samples or sample IDs. Is empty by default. If group assignments
#'             are given, the vector must have the same length as the number of
#'             samples. If sample IDs are provided the assignment will result
#'             in a binary group assignemt.
#' @param groupColSet A color set from RColorBrewer or a manual vector of
#'             colors, which length must match the number of categories
#'             from groups.
#' @param pch Integer or character to be used for plotting the points
#' @param col Set color for the points. If set, it must be a character vector
#'             of length 2. (1. normal point; 2. outlier point)
#' @param basePlot if \code{TRUE}, use the R base plot version, else use the
#'             plotly framework, which is the default
#' @param legendPos Set legendpos, by default topleft.
#' @param conf.alpha If set, a confidence interval is plotted, defaults to 0.05
#' @param samplePoints Sample points for Q-Q plot, defaults to max 30k points
#' @param xlim,ylim The x/y limits for the plot or NULL to use
#'             the full data range
#' @param ymax If set, ymax is the upper bound for the plot range on the y axis.
#' @param log If TRUE, the default, counts are plotted in log10.
#' @param rowCentered If TRUE, the counts are row-wise (gene-wise) centered
#' @param rowGroups,colGroups A vector of co-factors (colnames of colData)
#'             for color coding the rows. It also accepts a data.frame of
#'             dim = (#samples, #groups). Must have more than 2 groups.
#' @param rowColSet,colColSet A color set from RColorBrewer/colorRampPalette
#' @param nRowCluster,nColCluster Number of clusters to show in the row and
#'             column dendrograms. If this argument is set the resulting
#'             cluster assignments are added to the OutriderDataSet.
#' @param show_names character string indicating whether to show 'none', 'row',
#'             'col', or 'both' names on the heatmap axes.
#' @param bcvQuantile quantile for choosing the cutoff for the biological 
#'             coefficient of variation (BCV)
#' @param nGenes upper limit of number of genes (defaults to 500). Subsets the
#'             top n genes based on the BCV.
#' @param nBreaks number of breaks for the heatmap color scheme. Default to 50.
#' @param bins Number of bins used in the histogram. Defaults to 100.
#' @param yadjust Option to adjust position of Median and 90 percentile labels.
#' @param ylab The y axis label
#' @param labCex The label cex parameter
#' @param labLine Option to move axis labels
#'
#### Additional ... parameter
#' @param ... Additional parameters passed to plot() or plot_ly() if not stated
#'             otherwise in the details for each plot function
#'
#' @details
#'
#' \code{plotAberrantPerSample}: The number of aberrant events per sample are
#' plotted sorted by rank. The ... parameters are passed on to the
#' \code{\link{aberrant}} function.
#'
#' \code{plotVolcano}: the volcano plot is sample-centric. It plots for a given
#' sample the negative log10 nominal P-values against the Z-scores for all
#' genes.
#'
#' \code{plotExpressionRank}: This function plots for a given gene the
#' expression level against the expression rank for all samples. This can
#' be used with normalized and unnormalized expression values.
#'
#' \code{plotQQ}: the quantile-quantile plot for a given gene or if
#' \code{global} is set to \code{TRUE} over the full data set. Here the
#' observed P-values are plotted against the expected ones in the negative
#' log10 space.
#'
#' \code{plotExpectedVsObservedCounts}: A scatter plot of the observed counts
#' against the predicted expression for a given gene. 
#'  
#' \code{plotCountCorHeatmap}: The correlation heatmap of the count data
#' of the full data set. Default the values are log transformed and
#' row centered. This function returns an OutriderDataSet with annotated
#' clusters if requested. The ... arguments are passed to the
#' \code{\link[pheatmap]{pheatmap}} function.
#'
#' \code{plotCountGeneSampleHeatmap}: A gene x sample heatmap of the raw or
#' normalized counts. By default they are log transformed and row centered.
#' Only the top 500 viable genes based on the BCV (biological coefficient 
#' of variation) is used by default. 
#' 
#' \code{plotSizeFactors}: The sizefactor distribution within the dataset. 
#' 
#' \code{plotFPKM}: The distribution of FPKM values. If the OutriderDataSet
#' object contains the \code{passedFilter} column, it will plot both FPKM
#' distributions for the expressed genes and for the filtered genes.
#'
#' \code{plotExpressedGenes}: A summary statistic plot on the number of genes
#' expressed within this dataset. It plots the sample rank (based on the 
#' number of expressed genes) against the accumulated statistics up to the 
#' given sample. 
#' 
#' \code{plotDispEsts}: Plots the dispersion of the OutriderDataSet
#' model against the normalized mean count. If autoCorrect is used it will also
#' estimate the dispersion without normalization for comparison.
#'
#' \code{plotPowerAnalysis}: The power analysis plot should give the user a
#' ruff estimate of the events one can be detected with OUTRIDER. Based on
#' the dispersion of the provided OUTRIDER data set the theoretical P-value
#' over the mean expression is plotted. This is done for different expression
#' levels. The curves are smooths to make the reading of the plot easier.
#' 
#' \code{plotEncDimSearch}: Visualization of the hyperparameter optimization. 
#' It plots the encoding dimension against the achieved loss (area under the 
#' precision-recall curve). From this plot the optimum should be choosen for
#' the \code{q} in fitting process. 
#' 
#' @return If base R graphics are used nothing is returned else the plotly or
#'             the gplot object is returned.
#'
#' @examples
#' ods <- makeExampleOutriderDataSet(dataset="Kremer")
#' implementation <- 'autoencoder'
#' \dontshow{
#'     # reduce the object size to speed up the calculations
#'     ods <- ods[1:400,1:80]
#'     implementation <- 'pca'
#' }
#'
#' mcols(ods)$basepairs <- 300 # assign pseudo gene length for filtering
#' ods <- filterExpression(ods)
#' ods <- OUTRIDER(ods, implementation=implementation)
#'
#' plotAberrantPerSample(ods)
#'
#' plotVolcano(ods, 49)
#' plotVolcano(ods, 'MUC1365', basePlot=TRUE)
#'
#' plotExpressionRank(ods, 35)
#' plotExpressionRank(ods, "NDUFS5", normalized=FALSE,
#'     log=FALSE, main="Over expression outlier", basePlot=TRUE)
#'
#' plotQQ(ods, 149)
#' plotQQ(ods, global=TRUE, outlierRatio=0.001)
#' 
#' plotExpectedVsObservedCounts(ods, 149)
#' plotExpectedVsObservedCounts(ods, "ATAD3C", basePlot=TRUE)
#' 
#' plotExpressedGenes(ods)
#' 
#' sex <- sample(c("female", "male"), dim(ods)[2], replace=TRUE)
#' colData(ods)$Sex <- sex
#' ods <- plotCountCorHeatmap(ods, nColCluster=4, normalized=FALSE)
#' ods <- plotCountCorHeatmap(ods, colGroup="Sex", colColSet="Set1")
#' table(colData(ods)$clusterNumber_4)
#' 
#' plotCountGeneSampleHeatmap(ods, normalized=FALSE)
#' plotCountGeneSampleHeatmap(ods, rowGroups="theta", 
#'         rowColSet=list(c("white", "darkgreen")))
#' 
#' plotSizeFactors(ods)
#' 
#' mcols(ods)$basepairs <- 1
#' mcols(ods)$passedFilter <- rowMeans(counts(ods)) > 10
#' plotFPKM(ods)
#'
#' plotDispEsts(ods, compareDisp=FALSE)
#'
#' plotPowerAnalysis(ods)
#'
#' \dontrun{
#' # for speed reasons we only search for 5 different dimensions
#' ods <- findEncodingDim(ods, params=c(3, 10, 20, 35, 50), 
#'         implementation=implementation)
#' plotEncDimSearch(ods)
#' }
#'
#' @rdname plotFunctions
#' @aliases plotFunctions plotVolcano plotQQ plotExpressionRank
#'             plotCountCorHeatmap plotAberrantPerSample plotFPKM
#'             plotDispEsts plotPowerAnalysis
#'

#' @rdname plotFunctions
#' @export
plotVolcano <- function(ods, sampleID, main, padjCutoff=0.05, zScoreCutoff=0,
                    pch=16, basePlot=FALSE, col=c("gray", "firebrick")){
    if(missing(sampleID)){
        stop("specify which sample should be plotted, sampleID = 'sample5'")
    }
    if(!all(c('padjust', 'zScore') %in% assayNames(ods))){
        stop('Calculate Z-scores and P-values first.')
    }
    if(is.logical(sampleID)){
        sampleID <- which(sampleID)
    }
    if(is.numeric(sampleID)){
        if(!(is.numeric(sampleID) && max(sampleID) <= ncol(ods))){
            stop('Sample index is out of bounds:',
                    paste(sampleID, collapse=", "))
        }
        sampleID <- colnames(ods)[sampleID]
    }
    if(!all(sampleID %in% colnames(ods))){
        stop("Sample ID is not in the data set.")
    }
    if(length(sampleID) > 1){
        ans <- lapply(sampleID, plotVolcano, ods=ods, padjCutoff=padjCutoff,
                zScoreCutoff=zScoreCutoff, basePlot=basePlot)
        return(ans)
    }
    if(missing(main)){
        main <- paste0("Volcano plot: ", sampleID)
    }

    if(is.null(rownames(ods))){
        rownames(ods) <- paste("feature", seq_len(nrow(ods)), sep="_")
    }

    dt <- data.table(
        GENE_ID   = rownames(ods),
        pValue    = pValue(ods)[,sampleID],
        padjust   = padj(ods)[,sampleID],
        zScore    = zScore(ods)[,sampleID],
        normCts   = counts(ods, normalized=TRUE)[,sampleID],
        medianCts = rowMedians(counts(ods, normalized=TRUE)),
        expRank   = apply(
                counts(ods, normalized=TRUE), 2, rank)[,sampleID],
        aberrant  = aberrant(ods, padjCutoff=padjCutoff,
                zScoreCutoff=zScoreCutoff)[,sampleID],
        color     = col[1])
    dt[aberrant == TRUE, color:=col[2]]

    # remove the NAs from the zScores for plotting
    dt[is.na(zScore),zScore:=0]

    p <- ggplot(dt, aes(zScore, -log10(pValue), color=color, text=paste0(
                "Gene ID: ", GENE_ID,
                "<br>Sample ID: ", sampleID,
                "<br>Median normcount: ", round(medianCts, 2),
                "<br>normcount: ", round(normCts, 2),
                "<br>expression rank: ", as.integer(expRank),
                "<br>nominal P-value: ", signif(pValue,3),
                "<br>adj. P-value: ", signif(padjust,3),
                "<br>Z-score: ", signif(zScore,2)))) + 
        geom_point() + 
        theme_bw() + 
        xlab("Z-score") + 
        ylab(expression(paste(-log[10], "(", italic(P), "-value)"))) + 
        ggtitle(main) + 
        scale_color_identity() + 
        theme(legend.position = 'none')
    
    if(isFALSE(basePlot)){
        p <- p + ylab(paste("-log<sub>10</sub>(<i>P</i>-value)"))
        return(ggplotly(p, tooltip="text"))        
    }
    p
}

#' @rdname plotFunctions
#' @export
plotQQ <- function(ods, geneID, main, global=FALSE, padjCutoff=0.05,
                    zScoreCutoff=0, samplePoints=TRUE, legendPos="topleft",
                    outlierRatio=0.001, conf.alpha=0.05,
                    pch=16, xlim=NULL, ylim=NULL, col=NULL){
    checkOutriderDataSet(ods)
    stopifnot(isScalarLogical(global))
    if(missing(geneID) & isFALSE(global)){
        stop('Please provide a geneID or set global to TRUE')
    }
    # Singel gene QQplot.
    if(isFALSE(global)){
        geneID <- getGeneIndex(geneID, ods)

        # Produce multiple qqplot if geneID is a vector.
        if(length(geneID)>1L){
            lapply(geneID, plotQQ, ods=ods, main=main, legendPos=legendPos,
                    col=col, global=FALSE)
            return(invisible())
        }
        #Plot QQplot for single gene.
        if(missing(main)){
            main <- paste0('Q-Q plot for gene: ', geneID)
        }
        if(is.null(col)){
            col <- c('black', 'firebrick')
        }
        pVal <- as.numeric(assay(ods[geneID,], 'pValue'))
        #plot all points with cex=1 for single gene.
        pointCex <- 1
        #data table with expected and observerd log10(pValues)
        aberrantEvent <- aberrant(ods[geneID,], padjCutoff=padjCutoff,
                zScoreCutoff=zScoreCutoff, by='sample')
        df <- data.table(obs= -log10(pVal), pch=pch, subset=FALSE,
                col=ifelse(aberrantEvent, col[2], col[1]))

    # global QQplot
    } else {
        if(missing(main)){
            main <- 'Global Q-Q plot'
        }
        if(is.null(col)){
            col <- c('#d95f02', '#1b9e77')
        }
        pVal <- as.numeric(assay(ods, 'pValue'))

        # Reducing Point size for global QQplot.
        pointCex <- .5

        #data table with expected and observerd log10(pValues)
        df <- data.table(obs= -log10(pVal), col=col[1], pch=pch, subset=FALSE)

        if(!is.null(outlierRatio)){
            odssub <- ods[,aberrant(ods, by='s', padjCutoff=padjCutoff,
                    zScoreCutoff=zScoreCutoff) < outlierRatio*length(ods)]
            if(ncol(odssub) > 0){
                pVal <- as.numeric(assay(odssub, 'pValue'))
                dfsub <- data.table(obs=-log10(pVal), col=col[2], pch=pch,
                        subset=TRUE)
                df <- rbind(df, dfsub)
            }
        }
    }

    # compute expected pValues.
    df <- df[order(subset, -obs)]

    # Correct p value if needed
    df[is.na(obs) | is.infinite(obs), obs:=1]
    minNonZeroP <- min(df[obs!=0, obs]) * 1e-2
    df[obs==0, obs:=minNonZeroP]

    df[,exp:=-log10(ppoints(.N)), by='subset']
    if(is.null(xlim)){
        xlim=range(df[,exp])
    }
    if(is.null(ylim)){
        ylim=range(df[,obs], na.rm=TRUE)
    }
    plot(NA, xlim=xlim, ylim=ylim, main=main,
            xlab=expression(
                    paste(-log[10], " (expected ", italic(P), "-value)")),
            ylab=expression(
                    paste(-log[10], " (observed ", italic(P), "-value)")))


    # confidence band
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    if(isTRUE(conf.alpha)){
        conf.alpha <- 0.05
    }
    if(is.numeric(conf.alpha)){
        exp <- df[subset==FALSE,exp]
        len <- length(exp)
        slen <- seq_len(len)
        getY <- function(x, exp){
            x1 <- exp[2]
            x2 <- exp[1]
            y1 <- -log10(x[2])
            y2 <- -log10(x[1])
            m <- (y2-y1)/(x2-x1)
            return(10^-(y1 + m*((x2+1)-x1)))
        }
        upper <- qbeta(    conf.alpha/2, slen, rev(slen))
        lower <- qbeta(1 - conf.alpha/2, slen, rev(slen))
        polygon(col="gray", border="gray", x=c(rev(exp), max(exp)+c(1,1), exp),
                y=-log10(c(
                        rev(upper), getY(upper, exp), getY(lower, exp), lower)))
    }

    #Add points to plot.
    if(isTRUE(samplePoints)){
        samplePoints <- c(5000, 30000)
    }
    if(is.numeric(samplePoints) & length(samplePoints) == 2){
        plotPoint <- df[,seq_len(.N) %in% c(seq_len(min(.N, samplePoints[1])),
                sample(seq_len(.N), size=min(.N, samplePoints[2]))),
                by='subset'][,V1]
        df <- df[plotPoint]
    }
    df[,points(exp, obs, pch=pch, col=col, cex=pointCex)]

    # diagonal and grid
    abline(0,1,col="firebrick")
    grid()

    #Add legend
    if(isTRUE(global)){
        legenddt <- data.table(onlyFull=c(TRUE, FALSE, TRUE),
                text=c("Full data set", "Filtered data set",
                        paste0("CI (\u03B1 = ", signif(conf.alpha, 2), ")")),
                lty=1, lwd=6, col=c(col, "gray"))
        if(length(unique(df[,subset])) == 1){
            legenddt <- legenddt[onlyFull == TRUE]
        }
        legenddt[,legend(legendPos, text, lty=lty, lwd=lwd, col=col)]
    } else {
        if(is.numeric(conf.alpha)){
            legend(legendPos, lty=1, lwd=7, col="gray",
                    paste0("CI (\u03B1 = ", signif(conf.alpha, 2), ")"))
        }
    }

    return(invisible())
}

#' @rdname plotFunctions
#' @export
plotExpectedVsObservedCounts <- function(ods, geneID, main, basePlot=FALSE,
                    log=TRUE, groups=c(), groupColSet='Set1', ...){

    # check user input
    checkOutriderDataSet(ods)
    if(is.null(normalizationFactors(ods))) {
        stop('Normalized counts missing')
    }
    if(missing(geneID)){
        stop("Please Specify which gene should be plotted, geneID = 'geneA'")
    }
    geneID <- getGeneIndex(geneID, ods)
    if (missing(main)) {
        main <- paste("Predicted expression plot:", geneID)
    }

    ods <- ods[geneID]
    cnts <- data.table(
            feature_id = geneID,
            sampleID   = colnames(ods),
            observed   = as.vector(counts(ods)) + isTRUE(log),
            expected   = as.vector(normalizationFactors(ods)) + isTRUE(log),
            aberrant   = as.vector(aberrant(ods)))

    # group assignment
    if(length(groups) != ncol(ods)) {
        tmp_group <- logical(nrow(cnts))
        tmp_group[cnts$sampleID %in% groups] <- TRUE
        groups <- tmp_group
    }

    # rename NA groups
    groups[is.na(groups)] <- 'NA'
    cnts[, group := groups]

    g <- ggplot(cnts, aes(expected, observed, text=paste0(
                "Gene ID: ", feature_id, "<br>", 
                "Sample ID: ", sampleID, "<br>",
                "Raw count: ", observed, "<br>",
                "Expected count: ", round(expected, 2), "<br>"))) +
        theme_bw() +
        geom_abline(slope = 1, intercept = 0) +
        labs(title = main,
                x=paste('Expected counts', ifelse(isTRUE(log), '+ 1', '')),
                y=paste('Raw counts', ifelse(isTRUE(log), '+ 1', '')))

    if(isTRUE(log)) {
        g <- g + scale_x_log10() + scale_y_log10()
    }

    point_mapping <- aes()
    # distinguish whether groups are given or not
    if(uniqueN(cnts$group) > 1 ) { # more than 1 group given
        # set color of groups
        if (is.character(groupColSet)) {
            g  <- g + scale_color_brewer(palette = groupColSet)
        } else {
            g <- g + scale_color_manual(values = groupColSet)
        }
        point_mapping$colour <- substitute(group)
        g <- g + geom_point(point_mapping)
    } else { # no grouping given
        point_mapping$colour <- substitute(aberrant)
        g <- g + geom_point(point_mapping) +
            scale_color_manual(values = c("gray", "firebrick")) +
            theme(legend.position = 'none')
    }

    if (isTRUE(basePlot)) {
        return(g)
    }
    ggplotly(g, tooltip="text")
}


#' @rdname plotFunctions
#' @export
plotExpressionRank <- function(ods, geneID, main, padjCutoff=0.05,
                    zScoreCutoff=0, normalized=TRUE, basePlot=FALSE, log=TRUE,
                    col=c("gray", "firebrick"), groups=c(),
                    groupColSet='Accent'){
    # check user input
    checkOutriderDataSet(ods)
    if(isTRUE(normalized) & is.null(sizeFactors(ods))){
        stop("Please calculate the sizeFactors or normalization factors ",
                "before plotting normlized counts.")
    }
    if(missing(geneID)){
        stop("Please Specify which gene should be plotted, geneID = 'geneA'")
    }
    geneID <- getGeneIndex(geneID, ods)

    if(length(col) != 2){
        stop("Please provide two colors as a vector.")
    }
    if(length(groups) == 0){
        groups <- logical(ncol(ods))
    } else if(all(groups %in% colData(ods)[,'sampleID'])){
        group <- colData(ods)[,'sampleID'] %in% groups
    } else if(length(groups) == ncol(ods)){
        groups[is.na(groups)] <- 'NA'
    } else {
        stop("Please provide meaningfull input for 'groups'.")
    }

    # apply over each gene if vector
    if(length(geneID) > 1){
        ans <- lapply(geneID, plotExpressionRank, ods=ods,
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff,
                basePlot=basePlot, normalized=normalized)
        return(ans)
    }
    stopifnot(isScalarValue(geneID))

    # subset ods
    ods <- ods[geneID,]

    dt <- data.table(
        sampleID  = colnames(ods),
        counts    = as.vector(counts(ods, normalized=normalized)) + isTRUE(log),
        rawcounts = as.vector(counts(ods)),
        group     = groups)

    dt[, medianCts:= median(counts)]
    dt[, norm_rank:= rank(counts, ties.method = 'first')]
    if(all(c('pValue', 'padjust', 'zScore') %in% assayNames(ods))){
        dt[, pValue   := assay(ods, 'pValue')[1,]]
        dt[, padjust  := assay(ods, 'padjust')[1,]]
        dt[, zScore   := assay(ods, 'zScore')[1,]]
        dt[, aberrant := aberrant(ods, padjCutoff=padjCutoff,
                zScoreCutoff=zScoreCutoff)[1,]]
    } else {
        dt[,aberrant:=FALSE]
    }
    ylab <- paste0(ifelse(isTRUE(normalized), "Normalized", "Raw"),
            " counts", ifelse(isTRUE(log), " + 1", ""))
    if(missing(main)){
        main <- paste("Expression rank plot:", geneID)
    }

    # create ggplot object
    g <- ggplot(data=dt, aes(x = norm_rank, y = counts)) +
        labs(title = main, x = 'Sample rank', y = ylab) +
        theme_bw()

    if(isTRUE(log)){
        g <- g + scale_y_log10()
    }
    
    point_mapping <- aes(col = aberrant)
if(isFALSE(basePlot)){
        point_mapping$text <- dt[,paste0(
            "Gene ID: ", geneID,
            "<br>Sample ID: ", sampleID,
            ifelse(uniqueN(groups) == 1, "", paste0("<br>Group: ", group)),
            "<br>Median normcount: ", round(medianCts, digits = 1),
            "<br>rawcount: ", rawcounts,
            "<br>expression rank: ", round(norm_rank, digits = 1),
            ifelse(!'padjust' %in% colnames(dt), '',
                paste0("<br>adj. P-value: ", sprintf("%1.1E", padjust),
                    "<br>nominal P-value: ", sprintf("%1.1E", pValue),
                    "<br>Z-score: ", round(zScore, digits = 1)))
        )]
    }

    # switch color mode
    if(uniqueN(groups) > 1){
        point_mapping$fill <- substitute(group)
        ignoreAesTextWarning({
            g <- g + 
                geom_point(size=3, stroke=0.5, pch=21, mapping=point_mapping) +
                scale_fill_brewer(palette=groupColSet) +
                scale_color_manual(values= c('grey30', col[2]))
        })
    } else {
        ignoreAesTextWarning({
            g <- g + geom_point(size = 2, point_mapping) +
                scale_color_manual(values = col) +
                theme(legend.position = 'none')
        })
    }

    if(isTRUE(basePlot)){
        return(g)
    }
    return(ggplotly(g))
}

#' @rdname plotFunctions
#' @export
plotCountCorHeatmap <- function(ods, normalized=TRUE, rowCentered=TRUE,
                    rowGroups=NA, rowColSet=NA, colGroups=NA,
                    colColSet=NA, nRowCluster=4, nColCluster=4,
                    main="Count correlation heatmap", basePlot=TRUE, nBreaks=50,
                    show_names=c("none", "row", "col", "both"), ...) {
    
    checkDeprication(names2check=c("rowGroups"="rowCoFactor", 
            "colGroups"="colCoFactor", "nRow/nCol-Cluster"="nCluster", 
            "not used anymore"="dendrogram", "show_names"="names"), ...)
    
    if(!isTRUE(basePlot)){
        return(plotCountCorHeatmapPlotly(
            ods, normalized=normalized, rowCentered=rowCentered,
            rowGroups=rowGroups, rowColSet=rowColSet,
            colGroups=colGroups, colColSet=colColSet,
            nRowCluster=nRowCluster, nColCluster=nColCluster,
            main=main, ...))
    }

    # correlation
    fcMat <- as.matrix(log2(counts(ods, normalized=normalized) + 1))
    if(isTRUE(rowCentered)){
        fcMat <- fcMat - rowMeans(fcMat)
    }
    ctscor <- cor(fcMat, method="spearman")

    # extract annotation and set clustering if requested
    annotation_row <- getAnnoHeatmap(x=ods, matrix=ctscor, groups=rowGroups,
            nClust=nRowCluster, extractFun=colData)
    annotation_col <- getAnnoHeatmap(x=ods, matrix=ctscor, groups=colGroups,
            nClust=nColCluster, extractFun=colData)
    if(!is.null(annotation_row) & "nClust" %in% colnames(annotation_row)){
        colData(ods)[[paste0('clusterNumber_', nRowCluster)]] <- 
                annotation_row[,'nClust']
    }
    if(!is.null(annotation_col) & "nClust" %in% colnames(annotation_col)){
        colData(ods)[[paste0('clusterNumber_', nColCluster)]] <- 
            annotation_col[,'nClust']
    }
    
    show_names <- match.arg(show_names, several.ok=FALSE)

    plotHeatmap(
        ods, ctscor,
        annotation_row = annotation_row,
        annotation_col = annotation_col,
        rowColSet = rowColSet,
        colColSet = colColSet,
        main = main,
        show_names = show_names,
        breaks=seq(-1, 1, length.out = nBreaks),
        ...
    )
}

plotCountCorHeatmapPlotly <- function(x, normalized=TRUE, rowCentered=TRUE,
                    rowGroups=NA, rowColSet=NA, 
                    colGroups=NA, colColSet=NA, 
                    nCluster=4, main="Count correlation heatmap", ...){
    # correlation
    fcMat <- as.matrix(log2(counts(x, normalized=normalized) + 1))
    if(isTRUE(rowCentered)){
        fcMat <- fcMat - rowMeans(fcMat)
    }
    ctscor <- cor(fcMat, method="spearman")

    # dendogram and clusters
    clusters <- cutree(hclust(dist(ctscor)), nCluster)
    clustCol <- getXColors(factor(clusters), "Dark2")
    if(is.numeric(nCluster)){
        colData(x)$clusterNumber <- names(clustCol)
    }

    if(is.na(rowGroups) & is.na(colGroups)){
        rowColSet <- colColSet
    }
    if(!is.na(colGroups)){
        if(is.na(colColSet)){
            colColSet <- "Dark2"
        }
        colCols <- getXColors(colData(x)[,colGroups], colColSet)
    } else {
        colCols <- clustCol
    }
    if(!is.na(rowGroups)){
        if(is.na(rowColSet)){
            rowColSet <- "Dark2"
        }
        colRows <- getXColors(colData(x)[,rowGroups], rowColSet)
    } else {
        colRows <- clustCol
    }

    p <- heatmaply(ctscor, col_side_colors=colCols, 
            row_side_colors=colRows, ...)
    
    return(p)
}


##### plotCountGeneSampleHeatmap #####
#' @rdname plotFunctions
#' @export
plotCountGeneSampleHeatmap <- function(ods, normalized=TRUE, rowCentered=TRUE,
                    rowGroups=NA, rowColSet=NA, colGroups=NA,
                    colColSet=NA, nRowCluster=4, nColCluster=4,
                    main="Count Gene vs Sample Heatmap", bcvQuantile=0.9,
                    show_names=c("none", "col", "row", "both"),
                    nGenes=500, nBreaks=50, ...){

    if (is.null(rownames(ods))) {
        message("Missing row names, using row numbers.")
        rownames(ods) <- paste0('row', seq_row(ods))
    } else if (is.null(colnames(ods))) {
        message("Missing column names, using column numbers.")
        colnames(ods) <- paste0('col', seq_col(ods))
    }

    ## Take top nGenes variable genes
    bcv <- 1/sqrt(robustMethodOfMomentsOfTheta(counts(ods), 1e4, 1e-2))
    if("theta" %in% colnames(mcols(ods))){
        bcv <- 1/sqrt(theta(ods))
    }
    rowData(ods)$BCV <- bcv
    ods_sub <- ods[!is.na(bcv) & 
            bcv > quantile(bcv, probs=bcvQuantile, na.rm=TRUE),]

    # take the top n genes if specified
    if(!is.null(nGenes)){
        ods_sub <- ods_sub[rank(rowData(ods_sub)$BCV) <= nGenes,]
        main = paste0(main, "\n", nrow(ods_sub), " most variable genes")
    }

    # count matrix
    if(isTRUE(normalized)){
        # autoencoder normalised values
        fcMat <- as.matrix(log2(counts(ods_sub, normalized=TRUE) + 1))
    } else {
        # normalize using sizeFactors
        fcMat <- as.matrix(log2(counts(ods_sub, normalized=FALSE) + 1))
        fcMat <- t(t(fcMat)/estimateSizeFactorsForMatrix(fcMat))
    }

    # center matrix
    if(isTRUE(rowCentered)) {
        fcMat <- fcMat - rowMeans(fcMat)
    }

    # row annotations
    annotation_row <- getAnnoHeatmap(x=ods_sub, matrix=fcMat, 
            groups=rowGroups, nClust=nRowCluster, extractFun=rowData)
    annotation_col <- getAnnoHeatmap(x=ods_sub, matrix=t(fcMat), 
            groups=colGroups, nClust=nColCluster, extractFun=colData)

    show_names <- match.arg(show_names, several.ok=FALSE)

    plotHeatmap(
        ods, fcMat,
        annotation_row = annotation_row,
        annotation_col = annotation_col,
        rowColSet = rowColSet,
        colColSet = colColSet,
        main = main,
        show_names = show_names,
        breaks=seq(-10, 10, length.out = nBreaks),
        ...
    )
}

plotHeatmap <- function(ods, mtx, annotation_row=NULL, annotation_col=NULL,
                    rowColSet=NA, colColSet=NA, main = "Heatmap",
                    show_names = c("none", "col", "row", "both"),
                    annotation_colors=NA, breaks=NA, ...){

    if (is.null(rownames(ods))) {
        message("Missing row names, using row numbers.")
        rownames(ods) <- paste0('row', seq_row(ods))
    } else if (is.null(colnames(ods))) {
        message("Missing column names, using column numbers.")
        colnames(ods) <- paste0('col', seq_col(ods))
    }

    # collect colours
    row_colors <- getAnnoColors(colorSet=rowColSet, annotation=annotation_row)
    col_colors <- getAnnoColors(colorSet=colColSet, annotation=annotation_col)
    tmp_anno_colors <- c(row_colors, col_colors)
    
    # use user defined colors and override autogenerated once if needed
    if(!missing(annotation_colors) & !isScalarNA(annotation_colors)){
        tmp_anno_colors <- c(annotation_colors, tmp_anno_colors)
    }
    
    # remove duplicates and NULL colors
    annotation_colors <- tmp_anno_colors[!duplicated(names(tmp_anno_colors))]
    annotation_colors <- annotation_colors[
            !vapply(annotation_colors, is.null, logical(1))]
    if(length(annotation_colors) == 0){
        annotation_colors <- NA
    }

    # option to show names in heatmap
    show_names <- match.arg(show_names, several.ok=FALSE)

    print(pheatmap(
        mat = mtx,
        breaks = breaks,
        color = bluered(length(breaks)),
        main = main,
        border_color = NA,
        show_rownames = (show_names == 'row' | show_names == 'both'),
        show_colnames = (show_names == 'col' | show_names == 'both'),
        annotation_row = annotation_row,
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        labels_row = getNiceName(rownames(mtx), 12),
        labels_col = getNiceName(colnames(mtx), 12),
        ...
    ))

    return(invisible(ods))
}

#' @rdname plotFunctions
#' @export
plotAberrantPerSample <- function(ods, main, padjCutoff=0.05, zScoreCutoff=0,
                    outlierRatio=0.001, col=brewer.pal(3, 'Dark2')[c(1,2)],
                    yadjust=c(1.2, 1.2), labLine=c(3.5, 3), ymax=NULL,
                    ylab="#Aberrantly expressed genes", labCex=par()$cex, ...){

    if(missing(main)){
        main <- 'Aberrant Genes per Sample'
    }

    count_vector <- sort(aberrant(ods, by="sample", padjCutoff=padjCutoff,
            zScoreCutoff=zScoreCutoff, ...))
    ylim <- c(0.4, max(1, count_vector)*1.1)
    if(!is.null(ymax)){
        ylim[2] <- ymax
    }
    replace_zero_unknown <- 0.5
    ticks <- c(replace_zero_unknown, signif(10^seq(
            from=0, to=round(log10(max(1, count_vector))), by=1/3), 1))

    labels_for_ticks <- sub(replace_zero_unknown, '0', as.character(ticks))

    bp= barplot2(
            replace(count_vector, count_vector==0, replace_zero_unknown),
            log='y', ylim=ylim, names.arg='', xlab='', plot.grid=TRUE,
            grid.col='lightgray', ylab='', yaxt='n', border=NA, xpd=TRUE,
            col=col[(!count_vector <= max(1, length(ods)*outlierRatio)) + 1],
            main=main)

    n_names <- floor(length(count_vector)/20)
    xnames= seq_len(n_names*20)
    axis(side=1, at= c(0,bp[xnames,]), labels= c(0,xnames))
    axis(side=2, at=ticks, labels= labels_for_ticks, ylog=TRUE, las=2)

    # labels
    mtext('Sample rank', side=1, line=labLine[1], cex=labCex)
    mtext(ylab, side=2, line=labLine[2], cex=labCex)

    # legend and lines
    hlines = c(Median=ifelse(median(count_vector)==0, replace_zero_unknown,
            median(count_vector)) , Quantile90=quantile(
                    count_vector, 0.9, names=FALSE))
    color_hline= c('black','black')
    abline(h= hlines, col=color_hline)
    text(x=c(1,1), y= hlines*yadjust, col=color_hline, adj=0,
            labels=c('Median', expression(90^th ~ 'percentile')))

    box()
}


#' @rdname plotFunctions
#' @export
plotFPKM <- function(ods, bins=100){
    fpkm <- fpkm(ods)
    colors <- c("grey60","darkgreen")
    if(!'passedFilter' %in% colnames(mcols(ods))){
        message(paste0('To see the difference between the filtering ',
                'run first the filterExpression() function.'))
        passed <- rep(TRUE, nrow(ods))
        colors <- c("darkgreen")
    } else {
        passed <- mcols(ods)[['passedFilter']]
    }

    histdata <- data.table(melt(fpkm, value.name = 'fpkm'),
            'passedFilter'=rep(passed, dim(fpkm)[2]))

    if(any(histdata$fpkm == 0)){
        numZero <- sum(histdata$fpkm == 0)
        message(paste0(numZero, " sample-gene combinations are zero. This is ",
                signif(numZero/nrow(histdata)*100, 3), "% of the data"))
        histdata <- histdata[fpkm != 0]
    }

    p <- ggplot(histdata, aes(fpkm, fill=passedFilter), alpha=0.5) +
        geom_histogram(bins=bins) +
        scale_fill_manual(values=colors) + 
        theme_bw() +
        theme(legend.position = c(0.1, 0.9)) +
        scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
        labs(x='FPKM', y='Frequency') + 
        ggtitle("FPKM distribution")
    p
}


plotDispEsts.OUTRIDER <- function(object, compareDisp, xlim, ylim,
                    main="Dispersion estimates versus mean expression", ...){
    # validate input
    if(!'theta' %in% names(mcols(object))){
        stop('Fit OUTRIDER first by executing ods <- OUTRIDER(ods) ',
                'or ods <- fit(ods)')
    }
    if(missing(compareDisp)){
        compareDisp <- !is.null(E(object))
    }

    # disp from OUTRIDER
    odsVals <- getDispEstsData(object)
    legText <- c("OUTRIDER fit")
    legCol <- c("firebrick")
    nonAutoVals <- NULL

    if(isTRUE(compareDisp)){
        #fit OUTRIDER without AutoCorrect
        ods2 <- OutriderDataSet(countData = counts(object))
        ods2 <- estimateSizeFactors(ods2)
        ods2 <- fit(ods2)
        nonAutoVals <- getDispEstsData(ods2, odsVals$mu)
    }
    if(missing(xlim)){
        xlim=range(odsVals$xpred)
    }
    if(missing(ylim)){
        ylim=range(1/odsVals$disp)
    }
    # plot dispersion
    plot(NA, xlim=xlim, ylim=ylim, pch='.', log="yx", xlab='Mean expression',
            ylab='Dispersion', main=main)
    if(!is.null(nonAutoVals)){
        points(odsVals$mu, 1/nonAutoVals$disp, pch='.', col="black")
    }
    points(odsVals$mu, 1/odsVals$disp, pch='.', col='firebrick')

    # plot fit
    lines(odsVals$xpred, odsVals$ypred, lwd=2, col="firebrick")
    if(isTRUE(compareDisp)){
        lines(odsVals$xpred, nonAutoVals$ypred, lwd=2, col="black")
        legText <- c("before correction fit", "OUTRIDER fit")
        legCol <- c('black', 'firebrick')
    }

    legend("bottomleft", legText, col=legCol, pch=20, lty=1, lwd=3)
}


#' @rdname plotFunctions
#' @exportMethod plotDispEsts
setMethod("plotDispEsts", signature(object="OutriderDataSet"),
        plotDispEsts.OUTRIDER)


#' @rdname plotFunctions
#' @export
plotPowerAnalysis <- function(ods){
    dispfit <-getDispEstsData(ods)
    m <- 10^seq.int(0,4,length.out = 1E4)
    d <- 1/dispfit$fit(m)
    dt<-rbindlist(lapply(c(0,0.1,0.2,0.3,0.5, 2,5,10), function(frac)
        data.table(mean=m, disp=d, frac=frac,
                pVal=pmin(0.5, pnbinom(round(frac * m), mu = m, size=d),
                        1 - pnbinom(round(frac * m), mu = m, size=d) +
                                dnbinom(round(frac * m), mu = m, size=d))
        )))

    dt[,negLog10pVal:=-log10(pVal)]
    dt[,Fraction:=as.factor(frac)]
    dt[,ExprType:= ifelse(frac<1, 'Downregulation', 'Overexpression')]
    ggplot(dt, aes(mean, negLog10pVal, col=Fraction, linetype=ExprType)) +
        geom_smooth(method=lm, formula = y ~ bs(x, 10), se = FALSE) +
        scale_x_log10(breaks=c(1,5,10,50,100,500,1000,5000,10000)) +
        labs(x="Mean", y='-log10(P-value)',color='Expression level',
                linetype='Expression type') + ylim(0,15)
}

#' @rdname plotFunctions
#' @export
plotEncDimSearch <- function(ods){
    if(is(ods, 'OutriderDataSet')){
        if(!'encDimTable' %in% colnames(metadata(ods)) &
                !is(metadata(ods)$encDimTable, 'data.table')){
            stop('Please run first the findEncodingDim before ',
                    'plotting the results of it.')
        }
        dt <- metadata(ods)$encDimTable
        q <- getBestQ(ods)
    } else {
        dt <- ods
        dt <- dt[, opt:=
                encodingDimension[which.max(evaluationLoss)[1]], by=zScore]
        q <- dt[opt == encodingDimension, opt]
    }

    if(!is.data.table(dt)){
        stop('Please provide the encDimTable from the OutriderDataSet object.')
    }
    if(!'zScore' %in% colnames(dt)){
        dt[,zScore:='Optimum']
        dt[,opt:=q]
    }
    dtPlot <- dt[,.(enc=encodingDimension, z=as.character(zScore),
            eva=evaluationLoss, opt)]
    ggplot(dtPlot, aes(enc, eva, col=z)) +
        geom_point() +
        scale_x_log10() +
        geom_smooth(method='loess') +
        ggtitle('Search for best encoding dimension') +
        geom_vline(data=dtPlot[opt == enc], show.legend=TRUE,
                aes(xintercept=enc, col=z, linetype='Optimum')) +
        geom_text(data=dtPlot[opt == enc], aes(y=0.0, enc-0.5, label=enc)) +
        labs(x='Encoding dimensions',
                y='Evaluation loss', col='Z score', linetype='Best Q') +
        scale_linetype_manual(values="dotted")
}

#' @rdname plotFunctions
#' @export
plotExpressedGenes <- function(ods, main='Statistics of expressed genes'){
    # labels and col names
    exp_genes_cols <- c(
        'sampleID'                         = 'sampleID',
        'Expressed\ngenes'                 = 'expressedGenes',
        'Union of\nexpressed genes'        = 'unionExpressedGenes',
        'Intersection of\nexpressed genes' = 'intersectionExpressedGenes',
        'Genes passed\nfiltering'          = 'passedFilterGenes',
        'Rank'                             = 'expressedGenesRank')

    # validate input
    if(!all(exp_genes_cols %in% names(colData(ods)))){
        stop('Compute expressed genes first by executing \n\tods <- ',
                'filterExpression(ods, addExpressedGenes=TRUE)')
    }

    dt <- as.data.table(colData(ods)[, exp_genes_cols])
    colnames(dt) <- names(exp_genes_cols)
    melt_dt <- melt(dt, id.vars = c('sampleID', 'Rank'))

    ggplot(melt_dt, aes(Rank, value)) +
        geom_point(aes(col = variable), size=1) +
        geom_line(aes(col = variable)) +
        theme_bw(base_size = 14) +
        ylim(0, NA) +
        theme(legend.position = 'top', legend.title = element_blank()) +
        labs(y = 'Number of genes', x = 'Sample rank', title = main) +
        scale_color_brewer(palette = 'Set1')
}

#' @rdname plotFunctions
#' @export
plotSizeFactors <- function(ods, basePlot=TRUE){
    checkSizeFactors(ods)
    plotdt <- data.table(sizeFactor=sizeFactors(ods),
            sampleID=colnames(ods))[order(sizeFactor)][,.(
                    sizeFactor, sampleID, rank=.I)]
    g <- ggplot(plotdt, aes(rank, sizeFactor, text=paste0(
        "SampleID: ", sampleID, "<br>",
        "Sizefactor: ", round(sizeFactor, 3), "<br>",
        "Rank: ", rank))) +
        geom_point() +
        theme_bw() +
        xlab("Sample rank") + 
        ggtitle("Size factor plot")
    
    if(isFALSE(basePlot)){
        g <- ggplotly(g, tooltip="text")
    }
    g
}


checkDeprication <- function(names2check, ...){
    lnames <- names(list(...))
    for(i in seq_along(names2check)){
        if(!names2check[i] %in% lnames){
            next
        }
        warning("We changed the API for the given parameter.\n",
                "  Please switch to the new parameter:\n",
                "\t'", names2check[i], "' --> '", names(names2check[i]), "'")
    }
}
