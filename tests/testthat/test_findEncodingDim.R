context("Testing Find Encoding Dimension: ")

# Test Find Encoding Dimension.
ods <- makeExampleOutriderDataSet(dataset = 'Kremer')
res <- findEncodingDim(ods)

test_that("Result has correct format", {
    expect_is(res, 'data.frame')
    expect_equal(names(res), c('params', 'evalLoss'))
})

# Test In silico outliers.
freq <- 1E-2
zScore=3
inj='both'
ods <- injectOutliers(ods, freq=freq, zScore=zScore, inj=inj)

test_that('In silico outliers',{
    expect_true(all(counts(ods)[
        assays(ods)[['trueCorruptions']]==0] == assay(ods, 'trueCounts')[
            assays(ods)[['trueCorruptions']]==0]))
    expect_true(any(counts(ods)[
        assays(ods)[['trueCorruptions']]!=0] != assay(ods, 'trueCounts')[
            assays(ods)[['trueCorruptions']]!=0]))
})
