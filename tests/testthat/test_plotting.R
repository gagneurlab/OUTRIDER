context("Testing plot functions: ")

test_that("plotting", {
    set.seed(42)
    ods <- makeExampleOutriderDataSet(n=40, m=40)
    
    mcols(ods)[['basepairs']] <- rnbinom(nrow(ods), mu=1e4, size=2)
    ods <- filterExpression(ods, filterGenes=FALSE, fpkmCutoff=1e3)
    ods <- OUTRIDER(ods, implementation='pca', q=2)
    sex <- sample(c("female", "male"), dim(ods)[2], replace=TRUE)
    colData(ods)$Sex <- sex
    
    
    expect_null(plotAberrantPerSample(ods))
    
    
    expect_equal(class(plotVolcano(ods, 30)), c("plotly", "htmlwidget"))
    expect_null(plotVolcano(ods, "sample_10", basePlot=TRUE))
    
    
    expect_equal(class(plotExpressionRank(ods, 23)), c("plotly", "htmlwidget"))
    expect_is(plotExpressionRank(ods, "feature_13", normalized=FALSE,
            log=FALSE, basePlot=TRUE), 'ggplot')
    
    
    expect_equal(class(plotExpectedVsObservedCounts(ods, 39)), 
            c("plotly", "htmlwidget"))
    expect_is(plotExpectedVsObservedCounts(ods, "feature_32", basePlot=TRUE),
            "ggplot")
    
    
    expect_is(plotExpressedGenes(ods), "ggplot")
    
    expect_null(plotQQ(ods, 1))
    expect_null(plotQQ(ods, "feature_20"))
    expect_null(plotQQ(ods, global=TRUE))
    expect_null(plotQQ(ods, global=TRUE, outlierRatio=0.001))
    expect_null(plotQQ(ods, global=TRUE, outlierRatio=NULL))
    
    
    expect_is(plotCountCorHeatmap(ods, colGroup="Sex", basePlot=FALSE),
            'plotly')
    ods <- plotCountCorHeatmap(ods)
    expect_invisible(plotCountCorHeatmap(ods, colGroup="Sex", normalized=FALSE))
    expect_invisible(plotCountCorHeatmap(ods, rowGroup="sampleID"))
    expect_true("clusterNumber_4" %in% colnames(colData(ods)))
    expect_is(ods, "OutriderDataSet")
    
    
    expect_invisible(
            plotCountGeneSampleHeatmap(ods, normalized=FALSE, bcvQuantile=0.5))
    expect_invisible(plotCountGeneSampleHeatmap(ods, rowGroups="theta", 
            rowColSet=list(c("white", "darkgreen"))))
    
    pass <- mcols(ods)$passedFilter
    # don't use mcols(ods)$XXX <- NULL since this is broken in R < 3.5
    mcols(ods) <- mcols(ods)[,!colnames(mcols(ods)) %in% "passedFilter"]
    expect_message(plotFPKM(ods), "To see the difference .*")
    mcols(ods)[,'passedFilter'] <- pass
    expect_is(plotFPKM(ods), 'ggplot')
    
    
    plotPowerAnalysis(ods)
    
})


