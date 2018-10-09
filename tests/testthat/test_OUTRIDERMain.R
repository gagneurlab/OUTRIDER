context("Testing OUTRIDER main function")

test_that("test main OUTRIDER function", {
    maxIter <- 6
    ods <- makeExampleOutriderDataSet(30, 30)
    ods <- estimateSizeFactors(ods)
    ods <- controlForConfounders(ods, iteration=maxIter)
    ods <- computePvalues(ods)
    ods <- computeZscores(ods)
    res <- results(ods, all=TRUE)
    expect_is(res, "data.table")
    
    ods2 <- OUTRIDER(ods, iteration=maxIter)
    res2 <- results(ods2, all=TRUE)
    
    expect_equal(ods, ods2)
    expect_equal(res, res2)
})
