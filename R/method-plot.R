

#' This function is defined because its used in a FraseR function, but should
#' not be used within this package!
#' @noRd
getPlotDistributionData <- function(gr=NULL, fds=NULL, type=NULL){
    stop("This method should not be called within the ods package!")
}

#'
#' Plot a QQ-plot for a given gene
#' 
#' @param ods OutriderDataSet
#' @param geneID Gene to be plotted
#' @param global Flag to plot a global QQ-plot, default FALSE
#' @param padj significance level to mark outliers
#' @param subset samples to subset for
#' @param main title for the plot
#' @param maxOutlier y-axis range for the plot
#' @param ... additional arguments for the internal plotQQ function.
#'             This arguments are currently used for development.
#' 
#' @return None
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(20, 400)
#' ods <- OUTRIDER(ods)
#' plotQQ(ods, 1)
#' 
#' @export
plotQQ <- function(ods, geneID=1:length(ods), global=FALSE, padj=0.05,
                    subset=1:ncol(ods), main=NULL, maxOutlier=NULL, ...){
    if(global == FALSE){
        if(!is.null(geneID)){
            ods <- ods[geneID,]
        }
        if(is.null(main)){
            main <- paste0('QQ-plot for gene: ', geneID)
        }
        sapply(1:length(ods), function(i){
            data <- list(pvalues=assays(ods)[['pValue']][i,])
            if(isScalarNumeric(padj, na.ok=FALSE)){
                data[['aberrant']] <- assays(ods)[['padjust']][i,] <= padj
            }
            plotQQplot.FraseR(data=data, maxOutlier=maxOutlier, main=main, ...)    
        })
    }
    else {
        if(is.null(main)){
            main <- 'Global QQ-plot'
        }
        pvalues <- assays(ods)[['pValue']][,subset]
        pvalues <- pvalues[
                sample(c(TRUE, FALSE), prob=c(0.01, 0.99)) | pvalues < 1E-4]
        data <- list(pvalues=pvalues)
        plotQQplot.FraseR(data=data, maxOutlier=maxOutlier, main=main, ...)
    }
    return(invisible())
}

#'
#' qqplot
#'
#' @noRd
plotQQplot.FraseR <- function(gr=NULL, fds=NULL, type=NULL, data=NULL, 
                    maxOutlier=2, conf.alpha=0.05, sample=FALSE, pch=16,
                    breakTies=FALSE, main="QQ-plot", legendPos="topleft"){
    if(isScalarLogical(conf.alpha)){
        conf.alpha <- ifelse(isTRUE(conf.alpha), 0.05, NA)
    }
    
    # get data
    if(is.null(data)){
        if(is.null(gr) | is.null(fds) | is.null(type)){
            stop(paste0("If data is not provided gr, fds and type", 
                    "needs to be passed on."))
        }
        #TODO data <- getPlotDistributionData(gr, fds, type)
        stop('Implement for FraseR!')
    }
    
    # points
    obs <- -log10(sort(data$pvalues))
    obs[is.na(obs)] <- 0
    if(length(obs) < 2){
        warning("There are no pvalues or all are NA!")
        return(FALSE)
    }
    if(breakTies){
        obs <- breakTies(obs, logBase=10, decreasing=TRUE)
    }
    exp <- -log10(ppoints(length(obs)))
    len <- length(exp)
    
    # limits for nice plotting
    maxPoint <- max(c(exp, obs))
    ylim <- range(0, maxPoint)
    if(isScalarNumeric(maxOutlier, na.ok=FALSE)){
        ylim <- range(0, min(exp[1]*maxOutlier, maxPoint))
    }
    
    # main plot area
    plot(NA, main=main, xlim=range(exp), ylim=ylim,
            xlab=expression(-log[10] ~  "(expected P-value)"),
            ylab=expression(-log[10] ~ "(observed P-value)"))
    
    
    # confidence band
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    if(is.numeric(conf.alpha)){
        getY <- function(x, exp){
            x1 <- exp[2]
            x2 <- exp[1]
            y1 <- -log10(x[2])
            y2 <- -log10(x[1])
            m <- (y2-y1)/(x2-x1)
            return(10^-(y1 + m*((x2+1)-x1)))
        }
        upper <- qbeta(conf.alpha/2,   1:len, rev(1:len))
        lower <- qbeta(1-conf.alpha/2, 1:len, rev(1:len))
        polygon(col="gray", border="gray", x=c(rev(exp), max(exp)+c(1,1), exp),
                y=-log10(c(
                        rev(upper), getY(upper, exp), getY(lower, exp), lower)))
    }
    
    # grid
    grid()
    
    plotPoint <- TRUE
    if(isTRUE(sample)){
        lo <- length(obs)
        plotPoint <- 1:lo %in% unique(c(1:min(lo, 100), sort(sample(1:lo,
                size=min(lo, 30000), prob=log10(1+lo:1)/sum(log10(1+lo:1))))))
    }
    
    # add the points
    #TODO- how to handle very extreme outliers?
    #TODO maybe we can set cex = 0.5 for the points in case of the global QQplot
    outOfRange <- obs > max(ylim)
    points(exp[plotPoint & !outOfRange], obs[plotPoint & !outOfRange], pch=pch, 
            cex=1)
    
    # diagonal and grid
    abline(0,1,col="firebrick")
    
    # plot outliers
    if(sum(outOfRange) > 0){
        points(exp[plotPoint & outOfRange], rep(max(ylim), sum(outOfRange)),
                pch=2, col='red')
    }
    
    if('aberrant' %in% names(data) && any(data$aberrant)){
        points(exp[1:sum(data$aberrant)], obs[1:sum(data$aberrant)], 
                pch=pch, col='firebrick')
    }
    
    if(is.numeric(conf.alpha)){
        legend(legendPos, paste0("CI (\u03B1 = ",
                signif(conf.alpha, 2), ")"), lty=1, lwd=7, col="gray")
    }
    return(invisible())
}

#'
#' breaks ties in a qq plot to get a better distributed p-value plot
#' @noRd
breakTies <- function(x, logBase=10, decreasing=TRUE){
    intervals <- sort(unique(c(0, x)))
    idxintervals <- findInterval(x, intervals)
    for(idx in as.integer(names(which(table(idxintervals) > 1)))){
        if(is.numeric(logBase)){
            minval <- logBase^-intervals[idx+1]
            maxval <- logBase^-intervals[idx]
            rand   <- runif(sum(idxintervals==idx), minval, maxval)
            rand   <- -log(rand, logBase)
        } else {
            minval <- intervals[idx]
            maxval <- intervals[idx+1]
            rand   <- runif(sum(idxintervals==idx), minval, maxval)
        }
        x[idxintervals==idx] <- rand
    }
    if(!is.na(decreasing)){
        x <- sort(x, decreasing=TRUE)
    }
}



#' 
#' Volcano plot
#' 
#' Volcano plot for a given sample over all genes.
#' 
#' @param ods OutriderDataSet
#' @param sampleID sample which should be plotted. 
#'        Can also be a vector of samples. 
#' @param padj padj cutoff
#' @param basePlot R base plot version of the plot.
#' @return None
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(10, 100)
#' ods <- OUTRIDER(ods)
#' plotVolcano(ods, 1)
#' 
#' @export 
plotVolcano <- function(ods, sampleID, padj=0.05, basePlot=FALSE){
    if(missing(sampleID)){
        stop("specify which sample should be plotted, sampleID = 'sample5'")
    }
    # TODO this will cause a bug if only an index is provided.
    # if(sampleID %in% colnames(ods)){
    #     stop('Sample not in data set')
    # }
    if(length(sampleID) > 1){
        sapply(sampleID, plotVolcano, ods=ods, padj=padj)
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
                    counts(ods, normalized=TRUE), 2, rank)[,sampleID])
    aberrant <- as.data.table(aberrant(ods[,sampleID]), padj)
    setnames(aberrant, 2, 'aberrant')
    dt <- cbind(dt, aberrant)
    
    # remove the NAs from the zScores for plotting
    dt[is.na(zScore),zScore:=0]
    
    colors <- c("gray", "red")
    if(all(dt[,aberrant==FALSE])){
        colors <- "gray"
    }
    if(basePlot == TRUE){
        dt[,plot(zScore, -log10(pValue), col=ifelse(aberrant, 'red', 'gray'), 
            pch=16, cex=.7, xlab='Z-score', ylab=expression(-log[10](P-value)))]
        grid()
        title(paste0("Volcano plot: ", sampleID))
    }else{
        plot_ly(
            data=dt,
            x=~zScore,
            y=~-log10(pValue),
            #y =~pValue,
            type="scatter",
            mode="markers",
            color=~aberrant,
            colors=colors,
            text=~paste0(
                "Gene ID: ", GENE_ID,
                "<br>Sample ID: ", sampleID,
                "<br>Median normcount: ", signif(medianCts, 2),
                "<br>normcount: ", signif(normCts, 2),
                "<br>expression rank: ", as.integer(expRank),
                "<br>nominal P-value: ", signif(pValue,3),
                "<br>adj. P-value: ", signif(padjust,3),
                "<br>Z-score: ", signif(zScore,2)
            )
        ) #%>% layout(yaxis = list(type = "log"))
    }
}



#' 
#' Expression rank plot
#' 
#' Plot expression over expression rank per gene.
#' The plot can be used before and after fitting. 
#' 
#' @param ods A OUTRIDER data set.
#' @param geneID Id of gene wich should be plotted. 
#' If multiple gene ids are supplied mulitple plots will be made.
#' @param padj padj cutoff for coloring significant genes
#' @param normalized if TRUE the normalized counts are used 
#'             otherwise the raw counts
#' @param basePlot R base plot version.
#' @return None
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(10, 100)
#' ods <- OUTRIDER(ods)
#' plotExpressionRank(ods, 1)
#' 
#' @export
plotExpressionRank <- function(ods, geneID, padj=0.05, normalized=TRUE, 
                    basePlot=FALSE){
    if(missing(geneID)){
        stop("Please Specify which gene should be plotted, geneID = 'geneA'")
    }
    #TODO this will cause a bug if only an index is provided.
    # if(geneID %in% rownames(ods)){
    #     stop('Gene not in data set')
    # }
    if(length(geneID) > 1){
        sapply(geneID, plotExpressionRank, ods=ods, padj=padj)
        return()
    }
    
    dt <- data.table(
        sampleID = colnames(ods),
        normcounts = as.integer(counts(ods, normalized=normalized)[geneID,]))
    
    dt[, medianCts:= median(normcounts)]
    dt[, norm_rank:= rank(normcounts, ties.method = 'first')]
    if('padjust' %in% assayNames(ods) & 'zScore' %in% assayNames(ods)){
        dt[, padjust:=assays(ods)[['padjust']][geneID,]]
        dt[, zScore:=assays(ods)[['zScore']][geneID,]]
        dt[, aberrant:=aberrant(ods, padj=padj)[geneID,]]
        colors <- ifelse(any(dt[,aberrant]), c("gray", "red"), "gray")
    } else {
        colors <- "gray"
    }

    if(basePlot == TRUE){
        dt[,plot(norm_rank, normcounts + 1, log = 'y', pch=16,
            col = ifelse(!aberrant, 'gray', 'red'),
            xlab = 'Sample rank', ylab = 'Normalized counts + 1', main=geneID)]
        grid()
    }else{
        plot_ly(
            data=dt,
            x=~norm_rank,
            y=~normcounts+1,
            type="scatter",
            mode="markers",
            color=~aberrant,
            colors=colors,
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
        ) %>% layout(xaxis=list(showline=TRUE, title='Sample Rank'),
                yaxis=list(type='log', dtick='D0', exponentformat='power',
                        title='Normalized counts + 1'))
    }
}


#' 
#' Correlation heatmap
#' 
#' Correlation heatmap of the count data of the given samples
#'
#' @param x An OutriderDataSet object
#' @param normalized if TRUE, the normalized counts are used, default.
#' @param rowCoFactor a vector of co-factors for color coding the rows
#' @param rowColSet A vector of colors or a color set from RColorBrewer
#' @param colCoFactor A vector of co-factors for color coding the columns
#' @param colColSet A vector of colors or a color set from RColorBrewer
#' @param nCluster An integer to be used for cutting the dendrogram into groups
#' @param main The title for the plot
#' @param annotateCluster if TRUE and nCluster is a integer the grouping 
#'             is saved in the returned OutriderDataSet
#' @param dendrogram character string indicating whether to draw 
#'             'none', 'row', 'column' or 'both' dendrograms.
#' @param names character string indicating whether to draw 
#'             'none', 'row', 'col', or 'both' names.
#' @param basePlot if TRUE the heatmap.2 function is used otherwise the
#'             interactive version from ggplotly is used.
#' @param ... additional arguments to the \code{heatmap.2} function 
#' @return OutriderDataSet object
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
plotCountCorHeatmap <- function(x, normalized=TRUE, rowCoFactor=NULL, 
                    rowColSet="Set1", colCoFactor=NULL, colColSet="Set2", 
                    nCluster=4, main="Count correlation heatmap", 
                    annotateCluster=TRUE, dendrogram='both', basePlot=TRUE,
                    names=c("both", "row", "col", "none"), ...){
    if(!isTRUE(basePlot)){
        return(plotCountCorHeatmapPlotly(x, normalized=normalized, 
                rowCoFactor=rowCoFactor, rowColSet=rowColSet, 
                colCoFactor=colCoFactor, colColSet=colColSet, 
                nCluster=nCluster, main=main, annotateCluster=annotateCluster, 
                dendrogram=dendrogram, ...))
    }
    
    colRows  <- NULL
    colCols  <- NULL
    clustCol <- NULL
    
    # correlation
    fcMat <- as.matrix(log2(counts(x, normalized=normalized) + 1))
    #TODO check what is wrong here - rowMeans seems to introduce an error.
    #ctscor <- cor(fcMat - rowMeans(fcMat), method="spearman")
    ctscor <- cor(fcMat, method="spearman")
    # nice names
    colnames(ctscor) <- substr(dimnames(ctscor)[[1]], 0, 12)
    rownames(ctscor) <- colnames(ctscor)
    
    # dendogram and clusters
    if(isScalarNumeric(nCluster) & nCluster > 0){
        clustCol <- getXColors(cutree(hclust(dist(ctscor)), nCluster))
    }
    if(annotateCluster & !is.null(clustCol)){
        colData(x)$clusterNumber <- names(clustCol)
    }
    
    # color bars
    if(is.null(rowCoFactor) & is.null(colCoFactor) & !is.null(colColSet)){
        rowColSet <- colColSet
    }
    if(!is.null(colCoFactor)){
        colCols <- getXColors(colData(x)[,colCoFactor], colColSet)
    } else {
        if(!is.null(nCluster) & !is.null(colColSet)){
            colCols <- clustCol
        }
    }
    if(!is.null(rowCoFactor)){
        colRows <- getXColors(colData(x)[,rowCoFactor], rowColSet)
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
    return(invisible(x))
}

#' @rdname plotCountCorHeatmap
plotCountCorHeatmapPlotly <- function(x, normalized=TRUE, rowCoFactor=NULL, 
                    rowColSet="Set1", colCoFactor=NULL, colColSet="Set2",
                    nCluster=4, main="Count correlation heatmap", 
                    annotateCluster=TRUE, dendrogram='both', ...){
    
    # correlation
    fcMat <- as.matrix(log2(counts(x, normalized=normalized) + 1))
    ctscor <- cor(fcMat - rowMeans(fcMat), method="spearman")
    
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
#' Plot aberrant events per sample
#' 
#' @rdname plotAberrantPerSample
#' @aliases plotAberrantPerSample plotAberrantPerSamplePlotly
#' 
#' @param ods OutriderDataSet
#' @param padj adjusted pvalue cutoff.
#' @param ... further arguments to \code{aberrant}
#' @return None
#' 
#' @examples
#' 
#' ods <- makeExampleOutriderDataSet()
#' ods <- OUTRIDER(ods)
#' 
#' plotAberrantPerSample(ods)
#' 
#' @export
plotAberrantPerSample <- function(ods, padj=0.05, ...){
    count_vector <- sort(aberrant(ods, by="sample", padj=padj,...))
    ylim = c(0.4, max(1, count_vector)*1.1)
    xlab_line= 3.5
    ylab_line= 3
    replace_zero_unknown = 0.5
    #ticks= c(replace_zero_unknown, 1,2,5, 10,20,50, 100,200,500),
    ticks= c(replace_zero_unknown, signif(10^seq(
            from=0, to=round(log10(max(1, count_vector))), by=1/3), 1))
    #yshift_hlines_text= c(0.1, 4)

    
    labels_for_ticks = sub(replace_zero_unknown, '0', as.character(ticks))
    
    bp= barplot2(
        replace(count_vector, count_vector==0, replace_zero_unknown),
        log='y', ylim= ylim, 
        names.arg='', xlab= '', plot.grid=TRUE, grid.col='lightgray',
        ylab= '', 
        yaxt = 'n', border=NA, xpd=TRUE,
        col = ifelse(count_vector < length(ods)*0.001,'#fdb462','grey')
    )
    n_names <- floor(length(count_vector)/20)
    xnames= c(1:n_names*20)
    axis(side=1, at= c(0,bp[xnames,]), labels= c(0,xnames))
    axis(side=2, at=ticks, labels= labels_for_ticks, ylog=TRUE, las=2)
    
    # labels
    mtext( 'Sample rank', side=1, line=xlab_line)
    mtext( '#Aberrantly expressed genes', side=2, line= ylab_line)
    
    # legend and lines
    hlines = c(Median=ifelse(median(count_vector)==0, replace_zero_unknown,
            median(count_vector)) , Quantile90=quantile(
                    count_vector,0.9, names=FALSE))
    color_hline= c('black','black')
    abline(h= hlines, col=color_hline)
    text(x=c(1,1), y= hlines*1.2, labels= c('Median', '90th percentile'),
            col=color_hline, adj=0)
    box()
}


#' 
#' Plot histogram of FPKM values 
#'
#' @param ods a OutriderDataSet object containing the reads and fpkm values
#' @return a ggplot object containing the FPKM plot
#' 
#' @export
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
plotFPKM <- function(ods){
    fpkm <- fpkm(ods)
    passed <- mcols(ods)[['passedFilter']] 
    
    histdata <- data.table(melt(fpkm, value.name = 'fpkm'),
            'passedFilter'=rep(passed, dim(fpkm)[2]))
    
    if(any(histdata$fpkm == 0)){
        numZero <- sum(histdata$fpkm == 0)
        message(paste(numZero, "sample-gene combinations are zero. This is",
                signif(numZero/nrow(histdata)*100, 3), "% of the data"))
        histdata <- histdata[fpkm != 0]
    }
    
    p <- ggplot(histdata, aes(fpkm, fill = passedFilter)) +
        geom_histogram(bins = 100) +
        scale_fill_manual(values = c("grey","darkgreen")) +
        theme(legend.position = c(0.1, 0.9)) +
        scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
        labs(x='FPKM', y='Frequency')
    p
}
    
#'
#' Dispersion estimation plot
#' 
#' Plotting the dispersion of the OutriderDataSet model against the normalized 
#' mean count.
#' 
#' @param object An OutriderDataSet
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
setMethod("plotDispEsts", signature(object="OutriderDataSet"), 
                    function(object, compareDisp=TRUE){
    # disp from OUTRIDER
    odsVals <- getDispEstsData(object)
    legText <- c("OUTRIDER fit")
    legCol <- c("firebrick")
    
    # plot dispersion
    heatscatter(odsVals$mu, 1/odsVals$disp, pch='.', log="yx", 
            xlab='Mean expression', ylab='Dispersion',
            main="Heatscatter of dispersion estimates")
    
    if(checkAutoCorrectExists() & isTRUE(compareDisp)){
        #fit OUTRIDER without AutoCorrect 
        ods2 <- OutriderDataSet(countData = counts(object))
        ods2 <- estimateSizeFactors(ods2)
        ods2 <- fit(ods2)
        nonAutoVals <- getDispEstsData(ods2, odsVals$mu)
        
        points(odsVals$mu, 1/nonAutoVals$disp, pch='.', col="firebrick")
    }
    points(odsVals$mu, 1/odsVals$disp, pch='.', col='black')
    
    # plot fit
    lines(odsVals$xpred, odsVals$ypred, lwd=2, col="black")
    if(checkAutoCorrectExists() & isTRUE(compareDisp)){
        lines(odsVals$xpred, nonAutoVals$ypred, lwd=2, col="firebrick")
        legText <- c("OUTRIDER fit", "autoCorrect fit")
        legCol <- c('firebrick', "black")
    }
    
    legend("bottomleft", legText, col=legCol, pch=20, lty=1, lwd=3)
    
})

getDispEstsData <- function(ods, mu=NULL){
    odsMu <- rowMeans(counts(ods, normalized=TRUE))
    if(is.null(mu)){
        mu <- odsMu
    }
    disp <- mcols(ods)$disp
    xidx <- c(0.001, 0.01, 0.5, seq(1, max(mu), 1), 
            max(mu) * c(1.1, 1.5, 2, 5, 10, 100, 1000))
    
    # fit linear model
    fit <- lm(1 / disp ~ 1 + I(mu^-1))
    pred <- predict(fit, list(mu=xidx))
    
    return(list(
        mu=odsMu,
        disp=disp,
        xpred=xidx,
        ypred=pred
    ))
}
