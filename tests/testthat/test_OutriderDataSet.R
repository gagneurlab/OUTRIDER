context("Testing OutriderDataSet-class")

test_that("create test data set", {
    expect_is(makeExampleOutriderDataSet(), "OutriderDataSet")
})

test_that("constructur for OutriderDataSet", {
    dds <- makeExampleDESeqDataSet()
    expect_is(OutriderDataSet(se=dds), "OutriderDataSet")
    expect_is(OutriderDataSet(countData=counts(dds)), "OutriderDataSet")
    expect_is(OutriderDataSet(countData=counts(dds), colData=colData(dds)), "OutriderDataSet")
    
    ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
    expect_is(ods, "OutriderDataSet")
    expect_null(show(ods))
    expect_true(validObject(ods))
    
})
