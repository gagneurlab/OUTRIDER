context("Testing Find Encoding Dimension: ")

test_that("Test findEncodingDim", {
    
    # Test Find Encoding Dimension.
    ods <- makeExampleOutriderDataSet(dataset = 'Kremer')
    ods <- ods[1:40, 1:40]
    ods <- filterExpression(ods, filterGenes=TRUE, minCounts=TRUE)
    countsbefore <- counts(ods)
    expect_equal(getBestQ(ods), NA_integer_)
    ods <- findEncodingDim(ods, params=3:5, iteration=2)
    
    expect_equal(countsbefore, counts(ods))
    expect_is(ods, 'OutriderDataSet')
    expect_equal(metadata(ods)[['optimalEncDim']], getBestQ(ods))
    expect_is(getBestQ(ods), "integer")
    expect_is(metadata(ods)[['encDimTable']], 'data.table')
    
})

test_that('In silico outliers',{
    ods <- makeExampleOutriderDataSet()
    freq <- 1E-2
    ods <- injectOutliers(ods, freq=freq, zScore=3, lnorm=TRUE, sdlog=3,
            inj='both')
    
    tCor <- assay(ods, 'trueCorruptions')
    tCts <- assay(ods, 'trueCounts')
    numExpect <- nrow(ods) * ncol(ods) * freq
    
    expect_equal(counts(ods)[tCor == 0], tCts[tCor == 0])
    expect_true(all(counts(ods)[tCor != 0] != tCts[tCor != 0]))
})
