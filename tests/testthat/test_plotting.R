context("Testing plot functions: ")

test_that("plotting", {
    set.seed(42)
    ods <- makeExampleOutriderDataSet(n=20, m=20)
    ods <- OUTRIDER(ods, implementation='pca')
    mcols(ods)[['basepairs']] <- rnbinom(nrow(ods), mu=1e4, size=25)
    
    expect_is(plotCountCorHeatmap(ods, basePlot=FALSE), 'plotly')
    ods <- plotCountCorHeatmap(ods)
    expect_true("clusterNumber" %in% colnames(colData(ods)))
    expect_is(ods, "OutriderDataSet")
    
    expect_null(plotQQ(ods, 1))
    expect_null(plotQQ(ods, global=TRUE))
    expect_null(plotQQ(ods, global=TRUE, outlierRatio=0.001))
    expect_null(plotQQ(ods, global=TRUE, outlierRatio=NULL))

    expect_null(plotAberrantPerSample(ods))
    
    expect_message(plotFPKM(ods), "To see the difference .*")
    expect_is(plotFPKM(ods), 'ggplot')
    ods <- filterExpression(ods, filterGenes=FALSE, fpkmCutoff=1e5)
    expect_is(plotFPKM(ods), 'ggplot')
    
    out <- plotDispEsts(ods)
    expect_is(out, 'list')
    expect_length(out, 2)
    
    expect_equal(class(plotVolcano(ods, 1)), c("plotly", "htmlwidget"))
    expect_null(plotVolcano(ods, 1, basePlot=TRUE))
    
    expect_equal(class(plotExpressionRank(ods, 1)), c("plotly", "htmlwidget"))
    expect_null(plotExpressionRank(ods, 1, basePlot=TRUE))
    
})
