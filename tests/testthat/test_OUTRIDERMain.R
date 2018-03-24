context("Testing OUTRIDER main function")

test_that("test main OUTRIDER function", {
    ods <- makeExampleOutriderDataSet(50)
    ods <- estimateSizeFactors(ods)
    ods <- fit(ods)
    ods <- computePvalues(ods)
    ods <- computeZscores(ods)
    res <- results(ods, all=TRUE)
    expect_is(res, "data.table")
    
    ods2 <- OUTRIDER(ods, autoControl = FALSE)
    res2    <- results(ods2, all=TRUE)
    
    expect_equal(ods, ods2)
    expect_equal(res, res2)
})
