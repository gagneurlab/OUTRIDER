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
#'   \item plotCountCorHeatmap()
#'   \item plotFPKM()
#'   \item plotDispEsts()
#'   \item plotPowerAnalysis()
#'   \item plotEncDimSearch()
#' }
#' 
#' For a detailed description of each plot function please see the details.
#' Most of the functions share the same parameters. 
#' 
#### Data specific parameters
#' @param ods An OutriderDataSet, which is used to extract the data for plotting
#' @param object An OutriderDataSet, see ods parameter
#' @param sampleID A sampleID, which should be plotted. Can also be a vector.
#' @param geneID A geneID, which should be plotted. Can also be a vector.
#' @param padjCutoff Significance level to mark outliers
#' @param zScoreCutoff Z-score cutoff to mark outliers
#' @param global Flag to plot a global Q-Q plot, default FALSE
#' @param outlierRatio The fraction to be used for the outlier sample filtering
#' @param normalized If TRUE, the normalized counts are used, the default,
#'             otherwise the raw counts
#' @param compareDisp If TRUE, the default, and if the autoCorrect normalization
#'             was used it computes the dispersion without autoCorrect and 
#'             plots it for comparison.
#### Graphical parameters
#' @param main Title for the plot, if missing a default title will be used.
#' @param pch Integer or character to be used for plotting the points
#' @param col Set color for the points. If set, it must be a character vector 
#'             of length 2. (1. normal point; 2. outlier point)
#' @param basePlot if TRUE, use the R base plot version, else use the plotly 
#'             framework, which is the default
#' @param legendPos Set legendpos, by default topleft.
#' @param conf.alpha If set, a confidence interval is plotted, defaults to 0.05
#' @param samplePoints Sample points for Q-Q plot, defaults to max 30k points
#' @param xlim The x limits for the plot or NULL to use the full data range
#' @param ylim The y limits for the plot or NULL to use the full data range
#' @param log If TRUE, the default, counts are plotted in log10.
#' @param rowCentered If TRUE, the counts are row-wise (gene-wise) centered
#' @param rowCoFactor A vector of co-factors for color coding the rows
#' @param rowColSet A vector of colors or a color set from RColorBrewer
#' @param colCoFactor A vector of co-factors for color coding the columns
#' @param colColSet A vector of colors or a color set from RColorBrewer
#' @param nCluster An integer to be used for cutting the dendrogram into groups.
#'             If this argument is set the resulting clusters are saved in the 
#'             returned OutriderDataSet.
#' @param dendrogram A character string indicating whether to draw 
#'             'none', 'row', 'column' or 'both' dendrograms.
#' @param names character string indicating whether to draw 
#'             'none', 'row', 'col', or 'both' names.
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
#' \code{plotCountCorHeatmap}: The correlation heatmap of the count data
#' of the full data set. Default the values are log transformed and 
#' row centered. This function returns an OutriderDataSet with annotated 
#' clusters if requested. The ... arguments are passed to the 
#' \code{\link[gplots]{heatmap.2}} function.
#' 
#' \code{plotFPKM}: The distribution of FPKM values. If the OutriderDataSet
#' object contains the \code{passedFilter} column, it will plot both FPKM
#' distributions for the expressed genes and for the filtered genes.
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
#' @return If base R graphics are used nothing is returned else the plotly or
#'             the gplot object is returned.
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(dataset="Kremer")
#' \dontshow{
#'     # reduce the object size to speed up the calculations
#'     ods <- ods[400:410,60:70]
#' }
#' ods <- filterExpression(ods, minCounts=TRUE)
#' ods <- OUTRIDER(ods)
#' 
#' plotAberrantPerSample(ods)
#' 
#' plotVolcano(ods, 1)
#' plotVolcano(ods, 'MUC1404', basePlot=TRUE)
#' 
#' plotExpressionRank(ods, 1)
#' plotExpressionRank(ods, "FAAH", normalized=FALSE, 
#'     log=FALSE, main="Over expression outlier", basePlot=TRUE)
#' 
#' plotQQ(ods, 1)
#' plotQQ(ods, global=TRUE, outlierRatio=0.001)
#' 
#' sex <- sample(c("female", "male"), dim(ods)[2], replace=TRUE)
#' colData(ods)$sex <- sex
#' ods <- plotCountCorHeatmap(ods, colCoFactor="sex", normalized=FALSE)
#' ods <- plotCountCorHeatmap(ods, nCluster=4)
#' head(colData(ods)$clusterNumber)
#' 
#' mcols(ods)$basepairs <- 1
#' mcols(ods)$passedFilter <- rowMeans(counts(ods)) > 10
#' plotFPKM(ods)
#' 
#' plotDispEsts(ods, compareDisp=FALSE)
#' 
#' plotPowerAnalysis(ods)
#' 
#' ods <- findEncodingDim(ods)
#' plotEncDimSearch(ods)
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
        if(isTRUE(basePlot)){
            return(invisible())
        }
        return(ans)
    }
    if(missing(main)){
        main <- paste0("Volcano plot: ", sampleID)
    }
    
    dt <- data.table(
        GENE_ID   = rownames(ods),
        pValue    = assays(ods)[['pValue']][,sampleID],
        padjust   = assays(ods)[['padjust']][,sampleID],
        zScore    = assays(ods)[['zScore']][,sampleID],
        normCts   = as.vector(counts(ods[,sampleID], normalized=TRUE)),
        medianCts = rowMedians(counts(ods, normalized=TRUE)),
        expRank   = apply(
                counts(ods, normalized=TRUE), 2, rank)[,sampleID],
        aberrant  = aberrant(ods, padjCutoff=padjCutoff,
                zScoreCutoff=zScoreCutoff)[,sampleID],
        color=col[1])
    dt[aberrant == TRUE, color:=col[2]]
    
    # remove the NAs from the zScores for plotting
    dt[is.na(zScore),zScore:=0]
    
    if(isTRUE(basePlot)){
        dt[,plot(zScore, -log10(pValue), col=color, pch=pch, cex=.7, 
                main=main, xlab='Z-score', ylab=expression(
                        paste(-log[10], "(", italic(P), "-value)")))]
        grid(equilogs=FALSE)
        return(invisible())
    }
    plot_ly(
        data=dt,
        x=~zScore,
        y=~-log10(pValue),
        #y =~pValue,
        type="scatter",
        mode="markers",
        marker = list(
            color=~color
        ),
        text=~paste0(
            "Gene ID: ", GENE_ID,
            "<br>Sample ID: ", sampleID,
            "<br>Median normcount: ", round(medianCts, 2),
            "<br>normcount: ", round(normCts, 2),
            "<br>expression rank: ", as.integer(expRank),
            "<br>nominal P-value: ", signif(pValue,3),
            "<br>adj. P-value: ", signif(padjust,3),
            "<br>Z-score: ", signif(zScore,2)
        )
    ) %>% layout(title=main, 
            yaxis=list(title="-log<sub>10</sub>(<i>P</i>-value)"))
}

#' @rdname plotFunctions
#' @export
plotQQ <- function(ods, geneID, main, global=FALSE, padjCutoff=0.05, 
                zScoreCutoff=0, samplePoints=TRUE, 
                legendPos="topleft", outlierRatio=0.001, conf.alpha=0.05, 
                pch=16, xlim=NULL, ylim=NULL, col=NULL){
    if(!is(ods, 'OutriderDataSet')){
        stop('Please provide an OutriderDataSet')
    }
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
            
                dfsub <- data.table(obs= -log10(pVal), col=col[2], pch=pch,
                        subset=TRUE)
                df <- rbind(df, dfsub)
            }
        }
    }
    
    # compute expected pValues.
    df <- df[order(subset, -obs)]
    df[,exp:=-log10(ppoints(.N)), by='subset']
    if(is.null(xlim)){
        xlim=range(df[,exp])
    }
    if(is.null(ylim)){
        ylim=range(df[,obs])
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
plotExpressionRank <- function(ods, geneID, main, padjCutoff=0.05, zScoreCutoff=0, 
                    normalized=TRUE, basePlot=FALSE, log=TRUE,
                    col=c("gray", "firebrick"), ...){
    # check user input
    if(!is(ods, "OutriderDataSet")){
        stop("Please provide an OutriderDataSet")
    }
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
    
    # apply over each gene if vector
    if(length(geneID) > 1){
        ans <- lapply(geneID, plotExpressionRank, ods=ods, 
                padjCutoff=padjCutoff, zScoreCutoff=zScoreCutoff, 
                basePlot=basePlot, normalized=normalized)
        if(isTRUE(basePlot)){
            return(invisible())
        }
        return(ans)
    }
    
    # get data frame
    dt <- data.table(
        sampleID = colnames(ods),
        normcounts = as.integer(counts(ods, normalized=normalized)[geneID,]))
    dt[,normcounts:=normcounts + ifelse(isTRUE(log), 1, 0)]
    
    dt[, medianCts:= median(normcounts)]
    dt[, norm_rank:= rank(normcounts, ties.method = 'first')]
    dt[, color:=col[1]]
    if('padjust' %in% assayNames(ods) & 'zScore' %in% assayNames(ods)){
        dt[, padjust  := assays(ods)[['padjust']][geneID,]]
        dt[, zScore   := assays(ods)[['zScore']][geneID,]]
        dt[, aberrant := aberrant(ods, padjCutoff=padjCutoff,
                zScoreCutoff=zScoreCutoff)[geneID,]]
        dt[aberrant == TRUE, color:=col[2]]
    } else {
        dt[,aberrant:=FALSE]
    }
    ylab <- paste0(ifelse(isTRUE(normalized), "Normalized", "Raw"),
            " counts", ifelse(isTRUE(log), " + 1", ""))
    if(missing(main)){
        main <- geneID
    }
    
    # plot it
    if(isTRUE(basePlot)){
        dt[,plot(norm_rank, normcounts,
                log=ifelse(isTRUE(log), "y", ""), pch=16, col=color,
                main=main, xlab='Sample rank', ylab=ylab)]
        grid(equilogs=FALSE)
        return(invisible())
        
    }
    plot_ly(
        data=dt,
        x=~norm_rank,
        y=~normcounts,
        type="scatter",
        mode="markers",
        marker = list(
            color=~color
        ),
        text=~paste0(
            "Gene ID: ", geneID,
            "<br>Sample ID: ", sampleID,
            "<br>Median normcount: ", round(medianCts, digits = 1),
            "<br>normcount: ", round(normcounts, digits = 1),
            "<br>expression rank: ", round(norm_rank, digits = 1),
            if(any(names(assays(ods))== 'padjust')){
                paste0("<br>adj. P-value: ", sprintf("%1.1E", padjust))
            },
            if(any(names(assays(ods))== 'zScore')){
                paste0("<br>Z-score: ", round(zScore, digits = 1))
            }
        )
    ) %>% layout(title=main,
            xaxis=list(showline=TRUE, title='Sample Rank'),
            yaxis=list(type=ifelse(isTRUE(log), 'log', 'linear'), 
                    dtick=ifelse(isTRUE(log), "D1", 
                            signif(diff(range(dt[,normcounts]))/10, 1)), 
                    exponentformat='power',
                    title=ylab))
}


#' @rdname plotFunctions
#' @export
plotCountCorHeatmap <- function(ods, normalized=TRUE, rowCentered=TRUE, 
                    rowCoFactor=NULL, rowColSet="Set1", 
                    colCoFactor=NULL, colColSet="Set2", nCluster=4, 
                    main="Count correlation heatmap",
                    dendrogram='both', basePlot=TRUE,
                    names=c("both", "row", "col", "none"), ...){
    if(!isTRUE(basePlot)){
        return(plotCountCorHeatmapPlotly(ods, normalized=normalized, 
                rowCoFactor=rowCoFactor, rowColSet=rowColSet, 
                rowCentered=rowCentered,
                colCoFactor=colCoFactor, colColSet=colColSet, 
                nCluster=nCluster, main=main, dendrogram=dendrogram, ...))
    }
    
    colRows  <- NULL
    colCols  <- NULL
    clustCol <- NULL
    
    # correlation
    fcMat <- as.matrix(log2(counts(ods, normalized=normalized) + 1))
    if(isTRUE(rowCentered)){
        fcMat <- fcMat - rowMeans(fcMat)
    }
    ctscor <- cor(fcMat, method="spearman")
    # nice names
    colnames(ctscor) <- substr(dimnames(ctscor)[[1]], 0, 12)
    rownames(ctscor) <- colnames(ctscor)
    
    # dendogram and clusters
    if(isScalarNumeric(nCluster) & nCluster > 0){
        clustCol <- getXColors(cutree(hclust(dist(ctscor)), nCluster))
    }
    if(!is.null(clustCol)){
        colData(ods)$clusterNumber <- names(clustCol)
    }
    
    # color bars
    if(is.null(rowCoFactor) & is.null(colCoFactor) & !is.null(colColSet)){
        rowColSet <- colColSet
    }
    if(!is.null(colCoFactor)){
        colCols <- getXColors(colData(ods)[,colCoFactor], colColSet)
    } else {
        if(!is.null(nCluster) & !is.null(colColSet)){
            colCols <- clustCol
        }
    }
    if(!is.null(rowCoFactor)){
        colRows <- getXColors(colData(ods)[,rowCoFactor], rowColSet)
    } else {
        if(!is.null(nCluster)){
            colRows <- clustCol
        }
    }
    
    names <- match.arg(names, several.ok=FALSE)
    mpos <- c(1,1)
    if(names != "both"){
        if(names != "row"){
            colnames(ctscor) <- rep("", ncol(ctscor))
            mpos[1] <- 2.3
        }
        if(names != "col"){
            rownames(ctscor) <- rep("", nrow(ctscor))
            mpos[2] <- 2.3
        }
    }
    mpos <- c(7, 7) / mpos
    
    # assemble the arguments
    args <-list(x=ctscor, denscol='green', col=bluered(20), key.par=list(las=1),
            keysize=1, trace='none', key.ylab='', key.title='', main=main,
            margins=c(12, 11), key.xlab='Spearman corr.', dendrogram=dendrogram,
            add.expr=quote(mtext(c("Samples"),  c(1, 4),
                    tmpOUTRIDERMpos, cex=1.07, font=c(1))))
    if(!is.null(colRows)){
        args[["RowSideColors"]] <- colRows
    } 
    if(!is.null(colCols)){
        args[["ColSideColors"]] <- colCols
    }
    if(length(list(...)) > 0){
        args[names(list(...))] <- list(...)
    }
    env <- environment()
    env$heatmap.2 <- heatmap.2
    .GlobalEnv$tmpOUTRIDERMpos <- mpos
    
    # run it
    do.call("heatmap.2", args, quote = FALSE, envir = env)
    
    rm("tmpOUTRIDERMpos", envir=globalenv())
    return(invisible(ods))
}


plotCountCorHeatmapPlotly <- function(x, normalized=TRUE, rowCentered=TRUE,
                    rowCoFactor=NULL, 
                    rowColSet="Set1", colCoFactor=NULL, colColSet="Set2",
                    nCluster=4, main="Count correlation heatmap", 
                    dendrogram='both', ...){
    
    # correlation
    fcMat <- as.matrix(log2(counts(x, normalized=normalized) + 1))
    if(isTRUE(rowCentered)){
        fcMat <- fcMat - rowMeans(fcMat)
    }
    ctscor <- cor(fcMat, method="spearman")
    
    # nice names
    colnames(ctscor) <- substr(dimnames(ctscor)[[1]], 0, 12)
    rownames(ctscor) <- colnames(ctscor)
    
    # dendogram and clusters
    clustCol <- getXColors(cutree(hclust(dist(ctscor)), nCluster))
    if(is.numeric(nCluster)){
        colData(x)$clusterNumber <- names(clustCol)
    }
    
    if(is.null(rowCoFactor) & is.null(colCoFactor)){
        rowColSet <- colColSet
    }
    if(!is.null(colCoFactor)){
        colCols <- getXColors(colData(x)[,colCoFactor], colColSet)
    } else {
        colCols <- clustCol
    }
    if(!is.null(rowCoFactor)){
        colRows <- getXColors(colData(x)[,rowCoFactor], rowColSet)
    } else {
        colRows <- clustCol
    }
    
    clust <- ctscor %>% dist() %>% hclust()
    ord <- clust$order
    df <- ctscor[ord,]
    df <- df[,ord]
    p <- df %>%
        melt() %>% 
        ggplot(aes(Var1, Var2, fill = value)) + geom_tile()
    
    p <- plotly_build(p)
    p$data[[1]]$colorscale = list(c(0.9, "#ffe6cc"), c(1, "##ff0000"))
    
    p <- ggplotly(p)
    p
    return(p)
} 


#' @rdname plotFunctions
#' @export
plotAberrantPerSample <- function(ods, main, padjCutoff=0.05, zScoreCutoff=0,
                    outlierRatio=0.001,
                    col=brewer.pal(3, 'Dark2')[c(1,2)], yadjust=c(1.2, 1.2), 
                    labLine=c(3.5, 3), ylab="#Aberrantly expressed genes", 
                    labCex=par()$cex, ...){
    
    if(missing(main)){
        main <- 'Aberrant Genes per Sample'
    }
    
    count_vector <- sort(aberrant(ods, by="sample", padjCutoff=padjCutoff, 
            zScoreCutoff=zScoreCutoff, ...))
    ylim = c(0.4, max(1, count_vector)*1.1)
    replace_zero_unknown = 0.5
    ticks= c(replace_zero_unknown, signif(10^seq(
            from=0, to=round(log10(max(1, count_vector))), by=1/3), 1))
    
    labels_for_ticks = sub(replace_zero_unknown, '0', as.character(ticks))
    
    bp= barplot2(
        replace(count_vector, count_vector==0, replace_zero_unknown),
        log='y', ylim=ylim, names.arg='', xlab='', plot.grid=TRUE, 
        grid.col='lightgray', ylab='', yaxt='n', border=NA, xpd=TRUE,
        col=col[(!count_vector < length(ods)*outlierRatio) + 1], main=main)
    
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
                    count_vector,0.9, names=FALSE))
    color_hline= c('black','black')
    abline(h= hlines, col=color_hline)
    text(x=c(1,1), y= hlines*yadjust, col=color_hline, adj=0,
            labels=c('Median', expression(90^th ~ 'percentile')))
            
    box()
}


#' @rdname plotFunctions
#' @export
plotFPKM <- function(ods){
    fpkm <- fpkm(ods)
    if(!'passedFilter' %in% colnames(mcols(ods))){
        message(paste0('To see the difference between the filtering ', 
                'run first the filterExpression() function.'))
        passed <- rep(TRUE, nrow(ods))
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
    
    p <- ggplot(histdata, aes(fpkm, fill = passedFilter)) +
        geom_histogram(bins = 100) +
        scale_fill_manual(values = c("grey","darkgreen")) +
        theme(legend.position = c(0.1, 0.9)) +
        scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
        labs(x='FPKM', y='Frequency')
    p
}


plotDispEsts.OUTRIDER <- function(object, compareDisp, xlim, ylim, 
            main="Dispersion estimates versus mean expression", ...){
    # validate input                 
    if(!'disp' %in% names(mcols(object))){
        stop('Fit OUTRIDER first by executing ods <- OUTRIDER(ods) ',
                'or ods <- fit(ods)')
    } 
    if(missing(compareDisp)){
        compareDisp <- 'weights' %in% names(metadata(object))
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
                    dnbinom(round(frac * m), mu = m, size=d)
            )
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
    if(!'encDimTable' %in% colnames(metadata(ods)) & 
            !is(metadata(ods)$encDimTable, 'data.table')){
        stop('Please run first the findEncodingDim before ', 
                'plotting the results of it.')
    }
    
    ggplot(metadata(ods)$encDimTable, aes(encodingDimension, evaluationLoss)) +
        geom_point() + 
        scale_x_log10() + 
        geom_smooth(method='loess') +
        ggtitle('Search for best encoding dimension')
    
}
