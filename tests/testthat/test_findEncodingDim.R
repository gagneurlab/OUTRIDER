context("Testing Find Encoding Dimension: ")

test_that("Test findEncodingDim", {
    
    # Test Find Encoding Dimension.
    ods <- makeExampleOutriderDataSet(dataset = 'Kremer')
    expect_equal(getBestQ(ods), NA_integer_)
    ods <- findEncodingDim(ods)
    
    expect_is(ods, 'OutriderDataSet')
    expect_equal(metadata(ods)[['optimalEncDim']], 5)
    expect_equal(getBestQ(ods), 5)
    expect_is(metadata(ods)[['encDimTable']], 'data.table')
    
})

test_that('In silico outliers',{
    
    ods <- makeExampleOutriderDataSet(dataset = 'Kremer')
    freq <- 1E-2
    ods <- injectOutliers(ods, freq=freq, zScore=3, inj='both')
    
    tCor <- assay(ods, 'trueCorruptions')
    tCts <- assay(ods, 'trueCounts')
    numExpect <- nrow(ods) * ncol(ods) * freq
    
    expect_equal(counts(ods)[tCor == 0], tCts[tCor == 0])
    expect_true(all(counts(ods)[tCor != 0] != tCts[tCor != 0]))
    expect_gt(numExpect * 1.1, sum(tCor != 0))
    expect_lt(numExpect * 0.8, sum(tCor != 0))
    
})
