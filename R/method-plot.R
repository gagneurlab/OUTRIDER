#' 
#' Volcano plot
#' 
#' Volcano plot for a given sample over all genes.
#' 
#' @param ods An OutriderDataSet, which is used to extract the data for plotting
#' @param sampleID A sampleID, which should be plotted. Can also be a vector.
#' @param padjCutoff Significance level to mark outliers
#' @param zScoreCutoff Z-score cutoff to mark outliers
#' @param main Title for the plot, can be NULL, which is the default
#' @param pch Integer or character to be used for plotting the points
#' @param col Set color for the points. If set, it must be a character vector 
#'             of length 2. (1. normal point; 2. outlier point)
#' @param basePlot if TRUE, use the R base plot version, else use the plotly 
#'             framekwork, which is the default
#' @return The plotly object or NULL if if base R is used
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(1000, 100)
#' ods <- OUTRIDER(ods)
#' plotVolcano(ods, 1)
#' 
#' @export 
plotVolcano <- function(ods, sampleID, padjCutoff=0.05, zScoreCutoff=0, pch=16,
                    main=NULL, basePlot=FALSE, col=c("gray", "firebrick")){
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
        sapply(sampleID, plotVolcano, ods=ods, padjCutoff=padjCutoff, 
                zScoreCutoff=zScoreCutoff, basePlot=basePlot)
        return()
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
                xlab='Z-score', ylab=expression(
                        paste(-log[10], "(", italic(P), "-value)")))]
        grid(equilogs=FALSE)
        title(ifelse(!is.null(main), main, paste0("Volcano plot: ", sampleID)))
    }else{
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
        ) %>% layout(yaxis = list(title = "-log<sub>10</sub>(<i>P</i>-value)"))
    }
}


#'
#' Q-Q plots 
#' 
#' Plot a Q-Q plot for a given gene or a Q-Q plot over the full data set.
#' 
#' @inheritParams plotVolcano
#' @param geneID A gene as character or index, which should be plotted.
#'         It can also be a vector. 
#' @param legendPos Set legendpos, by default topleft.
#' @param global Flag to plot a global Q-Q plot, default FALSE
#' @param conf.alpha If set, a confidence interval is plotted
#' @param outlierRatio The fraction to be used for the outlier sample filtering
#' @param samplePoints Sample points for Q-Q plot, defaults to max 30k points
#' @param xlim The x limits for the plot or NULL to use the full data range
#' @param ylim The y limits for the plot or NULL to use the full data range
#' 
#' @return None
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(500, 300)
#' ods <- OUTRIDER(ods)
#' plotQQ(ods, 1)
#' plotQQ(ods, global=TRUE, outlierRatio=0.001)
#' 
#' @export
plotQQ <- function(ods, geneID=NULL, global=FALSE, padjCutoff=0.05, 
                zScoreCutoff=0, main=NULL, samplePoints=TRUE, 
                legendPos="topleft", outlierRatio=0.001, conf.alpha=0.05, 
                pch=16, xlim=NULL, ylim=NULL, col=NULL){
    if(!is(ods, 'OutriderDataSet')){
        stop('Please provide an OutriderDataSet')
    }
    stopifnot(isScalarLogical(global))
    if(is.null(geneID) & isFALSE(global)){
        stop('Please provide a geneID or set global to TRUE')
    }
    # Singel gene QQplot.
    if(isFALSE(global)){
        if(is.null(geneID)){
            stop('Please provide a geneID')
        }
        if(is.logical(geneID)){
            geneID <- which(geneID)
        }
        if(is.numeric(geneID)){
            if(!(is.numeric(geneID) && max(geneID) <= nrow(ods))){
                stop(paste('Gene index is out of bounds:', 
                           paste(geneID, collapse=", ")))
            }
            geneID <- rownames(ods)[geneID]
        }
        if(!all(geneID %in% rownames(ods))){
            stop("Gene ID is not in the data set.")
        }
    
        # Produce multiple qqplot if geneID is a vector.
        if(length(geneID)>1L){
            sapply(geneID, plotQQ, ods=ods, main=main, legendPos=legendPos,
                    col=col, global=FALSE)
            return()
        }
        #Plot QQplot for single gene.
        if(is.null(main)){
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
        if(is.null(main)){
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
            pVal <- as.numeric(assay(odssub, 'pValue'))
            
            dfsub <- data.table(obs= -log10(pVal), col=col[2], pch=pch,
                    subset=TRUE)
            df <- rbind(df, dfsub)
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



#' 
#' Expression rank plot
#' 
#' Plot expression over expression rank per gene.
#' The plot can be used before and after fitting. 
#' 
#' @inheritParams plotVolcano
#' @inheritParams plotQQ
#' @param normalized If TRUE the normalized counts are used,
#'             otherwise the raw counts
#' @param log If TRUE, the default, counts are plotted in log10.
#' @param ... Additional parameters passed to plot() or plot_ly().
#' @return None or a plotly object
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(1000, 100)
#' ods <- OUTRIDER(ods)
#' plotExpressionRank(ods, 1)
#' plotExpressionRank(ods, 1, normalized=FALSE, log=FALSE, main="1. Gene")
#' 
#' @export
plotExpressionRank <- function(ods, geneID, padjCutoff=0.05, zScoreCutoff=0, 
                    normalized=TRUE, basePlot=FALSE, main=NULL, log=TRUE,
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
    if(is.logical(geneID)){
        geneID <- which(geneID)
    }
    if(is.numeric(geneID)){
        if(!(is.numeric(geneID) && max(geneID) <= nrow(ods))){
            stop(paste('Gene index is out of bounds:', 
                    paste(geneID, collapse=", ")))
        }
        geneID <- rownames(ods)[geneID]
    }
    if(!all(geneID %in% rownames(ods))){
        stop('The gene IDs are not within the data set.')
    }
    if(length(col) != 2){
        stop("Please provide two colors as a vector.")
    }
    
    # apply over each gene if vector
    if(length(geneID) > 1){
        sapply(geneID, plotExpressionRank, ods=ods, padjCutoff=padjCutoff, 
                zScoreCutoff=zScoreCutoff, basePlot=basePlot, 
                normalized=normalized)
        return()
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
    main <- ifelse(!is.null(main), main, geneID)
    
    # plot it
    if(isTRUE(basePlot)){
        dt[,plot(norm_rank, normcounts,
                log=ifelse(isTRUE(log), "y", ""), pch=16, col=color,
                main=main, xlab='Sample rank', ylab=ylab)]
        grid(equilogs=FALSE)
    } else {
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
}


#' 
#' Correlation heatmap
#' 
#' Correlation heatmap of the count data of the given samples
#'
#' @inheritParams plotVolcano
#' @param normalized If TRUE, the normalized counts are used, the default
#' @param rowCentered If TRUE, the counts are row-wise (gene-wise) centered
#' @param rowCoFactor A vector of co-factors for color coding the rows
#' @param rowColSet A vector of colors or a color set from RColorBrewer
#' @param colCoFactor A vector of co-factors for color coding the columns
#' @param colColSet A vector of colors or a color set from RColorBrewer
#' @param nCluster An integer to be used for cutting the dendrogram into groups
#' @param annotateCluster If TRUE and nCluster is an integer the grouping 
#'             is saved in the returned OutriderDataSet
#' @param dendrogram A character string indicating whether to draw 
#'             'none', 'row', 'column' or 'both' dendrograms.
#' @param names character string indicating whether to draw 
#'             'none', 'row', 'col', or 'both' names.
#' @param ... Additional arguments passed to the 
#'             \code{\link[gplots]{heatmap.2}} function
#' @return An OutriderDataSet object
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- plotCountCorHeatmap(ods, annotateCluster=TRUE)
#' 
#' sex <- sample(c("female", "male"), dim(ods)[2], replace=TRUE)
#' colData(ods)$sex <- sex 
#' plotCountCorHeatmap(ods, colCoFactor="sex")
#' 
#' @export
plotCountCorHeatmap <- function(ods, normalized=TRUE, rowCentered=TRUE, 
                    rowCoFactor=NULL, rowColSet="Set1", 
                    colCoFactor=NULL, colColSet="Set2", nCluster=4, 
                    main="Count correlation heatmap", annotateCluster=TRUE,
                    dendrogram='both', basePlot=TRUE,
                    names=c("both", "row", "col", "none"), ...){
    if(!isTRUE(basePlot)){
        return(plotCountCorHeatmapPlotly(ods, normalized=normalized, 
                rowCoFactor=rowCoFactor, rowColSet=rowColSet, 
                rowCentered=rowCentered,
                colCoFactor=colCoFactor, colColSet=colColSet, 
                nCluster=nCluster, main=main, annotateCluster=annotateCluster, 
                dendrogram=dendrogram, ...))
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
    if(annotateCluster & !is.null(clustCol)){
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
                    annotateCluster=TRUE, dendrogram='both', ...){
    
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
    if(annotateCluster){
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


#' 
#' Aberrant events per sample plot
#' 
#' Plot aberrant events per sample
#' 
#' @inheritParams plotVolcano
#' @inheritParams plotQQ
#' @param yadjust Option to adjust position of Median and 90 percentile labels. 
#' @param ylab The y axis label
#' @param labCex The label cex parameter
#' @param labLine Option to move axis labels
#' @param ... Further arguments to \code{\link{aberrant}}
#' @return None
#' 
#' @rdname plotAberrantPerSample
#' @aliases plotAberrantPerSample plotAberrantPerSamplePlotly
#' 
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet(500, 100)
#' ods <- OUTRIDER(ods)
#' 
#' plotAberrantPerSample(ods)
#' 
#' @export
plotAberrantPerSample <- function(ods, padjCutoff=0.05, zScoreCutoff=0,
                    main=NULL, outlierRatio=0.001,
                    col=brewer.pal(3, 'Dark2')[c(1,2)], yadjust=c(1.2, 1.2), 
                    labLine=c(3.5, 3), ylab="#Aberrantly expressed genes", 
                    labCex=par()$cex, ...){
    
    if(is.null(main)){
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


#' 
#' Distribution plot of FPKM values
#' 
#' Plots the distribution of FPKM values. It also can compare between
#' filtered and non filtered features. 
#'
#' @inheritParams plotVolcano
#' @return A ggplot object containing the FPKM plot
#'
#' @examples
#' ctsFile <- system.file('extdata', 'GTExSkinSmall.tsv', package='OUTRIDER')
#' ods <- OutriderDataSet(countData=read.table(ctsFile, check.names=FALSE))
#' annotation <- system.file('extdata', 'gencode.v19.genes.small.gtf.gz', 
#'         package='OUTRIDER')
#' 
#' # Filter, storing the fpkm values and not subsetting the ods object.
#' ods <- filterExpression(ods, annotation, filterGenes=FALSE, savefpkm=TRUE)
#' 
#' #Display the pre distribution of counts.
#' plotFPKM(ods)
#' 
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
    
#'
#' Dispersion estimation plot
#' 
#' Plotting the dispersion of the OutriderDataSet model against the normalized 
#' mean count.
#' 
#' @param object An OutriderDataSet object containing the fitted model
#' @param compareDisp If TRUE, the default, and if the autoCorrect normalization
#'             was used it computes the dispersion without autoCorrect and 
#'             plots it for comparison.
#' @return None
#' 
#' @examples 
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- fit(ods)
#' plotDispEsts(ods)
#' 
#' @exportMethod plotDispEsts
setMethod("plotDispEsts", signature(object="OutriderDataSet"), 
                    function(object, compareDisp=NULL){
    # validate input                 
    if(!'disp' %in% names(mcols(object))){
        stop('Fit OUTRIDER first by executing ods <- OUTRIDER(ods) ',
                'or ods <- fit(ods)')
    } 
    if(is.null(compareDisp)){
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
    
    # plot dispersion
    plot(NA, xlim=range(odsVals$mu), ylim=range(1/odsVals$disp),
            pch='.', log="yx", xlab='Mean expression', ylab='Dispersion',
            main="Dispersion estimates versus mean expression")
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
})



#' 
#' PowerAnalysis plot
#'
#' The power analysis plot should give the user a ruff estimate of the events
#' one can be detected with OUTRIDER. Based on the dispersion of the provided
#' OUTRIDER data set the theoretical P-value over the mean expression is
#' plotted. This is done for different expression levels. The curves are 
#' smooths to make the reading of the plot easier.
#'
#' @inheritParams plotVolcano
#' @return The ggplot object
#'
#' @examples
#' ods <- makeExampleOutriderDataSet()
#' ods <- estimateSizeFactors(ods)
#' ods <- fit(ods)
#' plotPowerAnalysis(ods)
#' 
#' @export
plotPowerAnalysis <- function(ods){
    dispfit <-getDispEstsData(ods)
    m <- 10^seq.int(0,4,length.out = 1E4)
    d <- 1/dispfit$fit(m)
    dt<-rbindlist(lapply(c(0,0.1,0.2,0.3,0.5, 2,5,10), function(frac) 
        data.table(mean=m, disp=d, frac=frac, 
            pVal=pmin(0.5, pnbinom(frac * m, mu = m, size=d),
                1 - pnbinom(frac * m, mu = m, size=d) + 
                    dnbinom(frac * m, mu = m, size=d)
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
