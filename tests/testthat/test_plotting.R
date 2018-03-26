context("Testing plot functions: ")

test_that("plotting", {
    ods <- makeExampleOutriderDataSet(n=50)
    ods <- OUTRIDER(ods)
    
    ods <- plotCountCorHeatmap(ods)
    expect_true("clusterNumber" %in% colnames(colData(ods)))
    expect_is(ods, "OutriderDataSet")
    
    expect_null(plotQQ(ods, 1))
    expect_equal(class(plotVolcano(ods, 1)), c("plotly", "htmlwidget"))
    expect_equal(class(plotExpressionRank(ods, 1)), c("plotly", "htmlwidget"))
    expect_null(plotVolcano(ods, 1, basePlot=TRUE))
    expect_null(plotExpressionRank(ods, 1, basePlot=TRUE))
    
})
