context("Testing OutriderDataSet-class")

test_that("create test data set", {
    ods <- makeExampleOutriderDataSet(40,20)
    expect_is(ods, "OutriderDataSet")
    expect_equal(dim(ods), c(40,20))
    expect_equal(dim(metadata(ods)[['trueOutliers']]), c(40,20))
})

test_that("constructur for OutriderDataSet", {
    se <- as(DESeq2::makeExampleDESeqDataSet(), "SummarizedExperiment")
    expect_is(OutriderDataSet(se=se), "OutriderDataSet")
    expect_is(OutriderDataSet(countData=assay(se, 'counts')),
            "OutriderDataSet")
    expect_is(OutriderDataSet(countData=assay(se, 'counts'), 
            colData=colData(se)), "OutriderDataSet")
    
    ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
    expect_is(ods, "OutriderDataSet")
    expect_null(show(ods))
    expect_true(validObject(ods))
    
})
