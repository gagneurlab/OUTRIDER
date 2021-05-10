context("Testing Outrider2DataSet-class")

test_that("create test data set", {
    ods <- makeExampleOutrider2DataSet(20,20)
    expect_is(ods, "Outrider2DataSet")
    expect_equal(dim(ods), c(20,20))
    expect_equal(dim(assay(ods, 'trueOutliers')), c(20,20))
    expect_equal(profile(ods), "other")
    expect_null(show(ods))
    expect_true(validObject(ods))
})

test_that("create test protrider data set", {
    ods <- makeExampleProtriderDataSet(20,20)
    expect_is(ods, "Outrider2DataSet")
    expect_equal(dim(ods), c(20,20))
    expect_equal(dim(assay(ods, 'trueOutliers')), c(20,20))
    expect_equal(profile(ods), "protrider")
})

test_that("constructur for Outrider2DataSet", {
    se <- as(DESeq2::makeExampleDESeqDataSet(), "SummarizedExperiment")
    rownames(se) <- paste0("gene", seq_row(se))
    colData(se)[["sampleID"]] <- colnames(se)
    expect_is(Outrider2DataSet(se=se), "Outrider2DataSet")
    expect_is(Outrider2DataSet(inputData=assay(se, 'counts'), 
                                profile="outrider"), "Outrider2DataSet")
    expect_is(Outrider2DataSet(inputData=assay(se, 'counts'), 
                                profile="outrider", colData=colData(se)), 
                            "Outrider2DataSet")
})
