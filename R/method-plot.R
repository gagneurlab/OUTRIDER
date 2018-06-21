
#'
#' Plot a QQ-plot for a given gene
#' 
#' @param ods OutriderDataSet
#' @param geneID Gene to be plotted
#' @param global Flag to plot a global QQ-plot, default FALSE
#' @param padj significance level to mark outliers
#' @param zScore Z-score cutoff to coloring outliers
#' @param subset samples to subset for
#' @param main title for the plot
#' @param maxOutlier y-axis range for the plot
#' @param filterOutliers If TRUE, in the global QQ plot the full data set will
#'             be compared against the for outlier sample filtered data set
#' @param conf.alpha if set, a confidence interval is plotted
#' @param outlierRatio The fraction to be used for the outlier sample filtering
#' @param sample sample points for QQplot, only used, when global==TRUE.
#' @param legendPos set legendpos, by default topleft.
#' @param col set color
#' @param ... additional arguments for the internal plotQQ function.
#'             This arguments are currently used for development.
#' 
#' @return None
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(20, 400)
#' ods <- OUTRIDER(ods)
#' plotQQ(ods, 1)
#' plotQQ(ods, global=TRUE, filterOutliers=TRUE)
#' 
#' @export
plotQQ <- function(ods, geneID=NULL, global=FALSE, padj=0.05, zScoreCutoff=3,
                main=NULL, sample=TRUE, legendPos="topleft", outlierRatio=0.001,
                conf.alpha=0.05, pch=16, col=ifelse(isTRUE(global),
                c('#1b9e77', '#d95f02'), c('black', 'firebrick')), ...){
    
#TODO- how to handle very extreme outliers?
#TODO maybe we can set cex = 0.5 for the points in case of the global QQplot
    
    stopifnot(isScalarLogical(global))
    if(length(col) == 2 & isTRUE(global)){
        col <- col[2:1]
    }
    if(is.null(main)){
        if(isTRUE(global)){
            
        } else {
            main <- paste0('Q-Q plot for gene: ', geneID)
        }
    }
    # Singel gene QQplot.
    if(isFALSE(global)){
        if(is.null(geneID)){
            stop('Please provide a geneID')
        }
        # Produce multiple qqplot if geneID is a vector.
        if(length(geneID)>1L){
            sapply(geneID, plotQQ, main=main, legendPos=legendPos, col=col[1],
                   global=FALSE)
            return()
        }
        #Plot QQplot for single gene.
        if(is.null(main)){
            main <- paste0('Q-Q plot for gene: ', geneID)
        }
        pVal <- as.numeric(assay(ods[geneID,], 'pValue'))
        #plot all points with cex=1 for single gene.
        plotPoint <- TRUE
        pointCex <- 1
        
        #TODO why does the col ifelse above not work?
        col <- c('black', 'firebrick')
        #data table with expected and observerd log10(pValues)
        df <- data.table(obs= -log10(pVal), col=ifelse(aberrant(ods[geneID,], 
            padj=padj,zScore=zScoreCutoff), col[2],col[1]),
            pch=pch, subset=FALSE, plotPoint=plotPoint)[order(-obs)]  
    }
    # global QQplot
    else {
        if(is.null(main)){
            main <- 'Global Q-Q plot'
        }
        pVal <- as.numeric(assay(ods, 'pValue'))
        plotPoint <- TRUE
        pointCex <- 1
        
        col <- c('#1b9e77', '#d95f02')
        #data table with expected and observerd log10(pValues)
        df <- data.table(obs= -log10(pVal), 
            col=col[1],
            pch=pch, subset=FALSE, 
            plotPoint=plotPoint)[order(-obs)]  
        
        if(!is.null(outlierRatio)){
            odssub <- ods[,aberrant(ods, by='s', padj=padj,
                zScore=zScoreCutoff) < outlierRatio*length(ods)]
            pVal <- as.numeric(assay(odssub, 'pValue'))
            
            dfsub <- data.table(obs= -log10(pVal), 
                col=col[2],
                pch=pch, subset=TRUE, 
                plotPoint=plotPoint)[order(-obs)] 
            df <- rbind(df, dfsub)
        }
        df <- df[order(-obs)]
        
        if(isTRUE(sample)){
            df[,plotPoint:= 1:.N %in% unique(sort(c(1:min(.N, 5000),
                sample(1:.N, size=min(.N, 30000)))))]
        }
            
        # Reducing Point size for global QQplot.
        pointCex <- .5
    }
    # compute expected pValues.
    df[,exp:= -log10(ppoints(.N)), by='subset']
    
    plot(NA, xlim=range(df[,exp]), ylim=range(df[,obs]), main=main,
         xlab=expression(paste(-log[10], " (expected ", italic(P), "-value)")),
         ylab=expression(paste(-log[10], " (observed ", italic(P), "-value)")))
    
    
    # confidence band
    # http://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
    if(is.numeric(conf.alpha)){
        exp <- df[subset==FALSE,exp]
        len <- length(exp)
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
    
    #Add legend
    if(isTRUE(global)){
        legend(legendPos, c("Full data set", "Filtered data set", paste0(
        "CI (\u03B1 = ", signif(conf.alpha, 2), ")")),
        lty=1, lwd=6, col=c(col, "gray"))
    } else {
        if(is.numeric(conf.alpha)){
            legend(legendPos, paste0("CI (\u03B1 = ",
            signif(conf.alpha, 2), ")"), lty=1, lwd=7, col="gray")
        }
    }
    
    #Add points to plot.
    points(df[,exp], df[,obs],  
         pch=df[,pch], col=df[,col])
    
    # diagonal and grid
    abline(0,1,col="firebrick")
    grid()
    return(invisible())
}




#' 
#' Volcano plot
#' 
#' Volcano plot for a given sample over all genes.
#' 
#' @param ods OutriderDataSet
#' @param sampleID sample which should be plotted. 
#'        Can also be a vector of samples. 
#' @param padjCut padj cutoff
#' @param zScoreCut Z-score cutoff
#' @param basePlot R base plot version of the plot.
#' @param col colors for the points in the form of c(non outliers, outliers)
#' @param main string passed to main (title) of the plot.
#' @return None
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(10, 100)
#' ods <- OUTRIDER(ods)
#' plotVolcano(ods, 1)
#' 
#' @export 
plotVolcano <- function(ods, sampleID, padjCut=0.05, zScoreCut=3, 
                main=NULL, basePlot=FALSE, col=c("gray", "firebrick")){
    if(missing(sampleID)){
        stop("specify which sample should be plotted, sampleID = 'sample5'")
    }
    if(is.logical(sampleID)){
        sampleID <- which(sampleID)
    }
    if(is.numeric(sampleID)){
        if(!(is.numeric(sampleID) && max(sampleID) <= ncol(ods))){
            stop(paste('Sample index is out of bounds:', 
                    paste(sampleID, collapse=", ")))
        }
        sampleID <- colnames(ods)[sampleID]
    }
    if(!all(sampleID %in% colnames(ods))){
        stop("Sample ID is not in the data set.")
    }
    if(length(sampleID) > 1){
        sapply(sampleID, plotVolcano, ods=ods, padjCut=padjCut, 
                zScoreCut=zScoreCut, basePlot=basePlot)
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
    aberrant <- as.data.table(aberrant(ods[,sampleID]), padjCut, zScoreCut)
    setnames(aberrant, 2, 'aberrant')
    dt <- cbind(dt, aberrant)
    
    # remove the NAs from the zScores for plotting
    dt[is.na(zScore),zScore:=0]
    
    
    if(all(dt[,aberrant==FALSE])){
        col <- col[1]
    }
    if(basePlot == TRUE){
        dt[,plot(zScore, -log10(pValue), col=ifelse(aberrant, col[2], col[1]), 
            pch=16, cex=.7, xlab='Z-score', 
            ylab=expression(paste(-log[10], "(", italic(P), "-value)")))]
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
            color=~aberrant,
            colors=col,
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
#' @param padjCut padj cutoff for coloring significant genes
#' @param zScoreCut Z-score cutoff for coloring significant genes
#' @param normalized if TRUE the normalized counts are used 
#'             otherwise the raw counts
#' @param basePlot R base plot version.
#' @param main string passed to main (title) of the plot.
#' @return None
#' 
#' @examples
#' ods <- makeExampleOutriderDataSet(10, 100)
#' ods <- OUTRIDER(ods)
#' plotExpressionRank(ods, 1)
#' 
#' @export
plotExpressionRank <- function(ods, geneID, padjCut=0.05, zScoreCut=3, 
                    normalized=TRUE, basePlot=FALSE, main=NULL){
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
    if(length(geneID) > 1){
        sapply(geneID, plotExpressionRank, ods=ods, padjCut=padjCut, 
                zScoreCut=zScoreCut, basePlot=basePlot, normalized=normalized)
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
        dt[, aberrant:=aberrant(ods, padj=padjCut, zScore=zScoreCut)[geneID,]]
        colors <- ifelse(any(dt[,aberrant]), c("gray", "firebrick"), "gray")
        
    } else {
        dt[,aberrant:=FALSE]
        colors <- "gray"
    }

    if(isTRUE(basePlot)){
        dt[,plot(norm_rank, normcounts + 1, log = 'y', pch=16,
                col=ifelse(dt[,aberrant], 'firebrick', 'grey'), 
                main=ifelse(!is.null(main), main, geneID), 
                xlab = 'Sample rank', ylab = 'Normalized counts + 1')]
        grid(equilogs=FALSE)
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
#' @param rowCentered if TRUE, the counts are rowwise (genewise) centered
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
plotCountCorHeatmap <- function(x, normalized=TRUE, rowCentered=TRUE, 
                    rowCoFactor=NULL, rowColSet="Set1", 
                    colCoFactor=NULL, colColSet="Set2", nCluster=4, 
                    main="Count correlation heatmap", annotateCluster=TRUE,
                    dendrogram='both', basePlot=TRUE,
                    names=c("both", "row", "col", "none"), ...){
    if(!isTRUE(basePlot)){
        return(plotCountCorHeatmapPlotly(x, normalized=normalized, 
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
    fcMat <- as.matrix(log2(counts(x, normalized=normalized) + 1))
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
#' Plot aberrant events per sample
#' 
#' @rdname plotAberrantPerSample
#' @aliases plotAberrantPerSample plotAberrantPerSamplePlotly
#' 
#' @param ods OutriderDataSet
#' @param padj adjusted pvalue cutoff.
#' @param outlierRatio The fraction to be used for the outlier sample filtering
#' @param ... further arguments to \code{aberrant}
#' @param main string passed to main (title) of the plot.
#' @param col set color.
#' @param yadjust option to adjust position of Median and 90percentile labels. 
#' @param ylab y-axis label
#' @param labCex label cex parameter 
#' @param labLine option to move axis labels
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
plotAberrantPerSample <- function(ods, padj=0.05, main=NULL, outlierRatio=0.001,
                    col=brewer.pal(3, 'Dark2')[c(1,2)], yadjust=c(1.2, 1.2), 
                    labLine=c(3.5, 3), ylab="#Aberrantly expressed genes", 
                    labCex=par()$cex, ...){
    
    if(is.null(main)){
        main <- 'Aberrant Genes per Sample'
    }
    
    count_vector <- sort(aberrant(ods, by="sample", padj=padj, ...))
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
    xnames= c(1:n_names*20)
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
        points(odsVals$mu, 1/nonAutoVals$disp, pch='.', col="firebrick")
    }
    points(odsVals$mu, 1/odsVals$disp, pch='.', col='black')
    
    # plot fit
    lines(odsVals$xpred, odsVals$ypred, lwd=2, col="black")
    if(isTRUE(compareDisp)){
        lines(odsVals$xpred, nonAutoVals$ypred, lwd=2, col="firebrick")
        legText <- c("before correction fit", "autoCorrect fit")
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
    xidx <- 10^(seq.int(log10(min(mu))-1, log10(max(mu))+0.1, length.out = 500))

    # fit DESeq2 parametric Disp Fit
    fit <- parametricDispersionFit(mu, 1/disp)
    pred <- fit(xidx)
    return(list(
        mu=mu,
        disp=disp,
        xpred=xidx,
        ypred=pred,
        fit=fit
    ))
}

#'
#' This function is not exported from DESeq2. Therefore we copied it over to
#' here. If DESeq2 will export this function we can use it instead.
#' 
#' TODO
#' 
#' @noRd
parametricDispersionFit <- function (means, disps){
    coefs <- c(0.1, 1)
    iter <- 0
    while (TRUE) {
        residuals <- disps/(coefs[1] + coefs[2]/means)
        good <- which((residuals > 1e-04) & (residuals < 15))
        suppressWarnings({
            fit <- glm(disps[good] ~ I(1/means[good]), 
                    family = Gamma(link = "identity"), 
                    start = coefs)
        })
        oldcoefs <- coefs
        coefs <- coefficients(fit)
        if (!all(coefs > 0)) 
            stop(simpleError("parametric dispersion fit failed"))
        if ((sum(log(coefs/oldcoefs)^2) < 1e-06) & fit$converged) 
            break
        iter <- iter + 1
        if (iter > 10) 
            stop(simpleError("dispersion fit did not converge"))
    }
    names(coefs) <- c("asymptDisp", "extraPois")
    ans <- function(q) coefs[1] + coefs[2]/q
    attr(ans, "coefficients") <- coefs
    ans
}

