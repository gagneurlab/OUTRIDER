context("Testing Find Encoding Dimension: ")

test_that("Test findEncodingDim OUTRIDER", {
    
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

# helper function to skip tests if we don't have the 'py_outrider' module
skip_if_no_py_outrider <- function() {
    have_py_outrider <- reticulate::py_module_available("py_outrider")
    if (!have_py_outrider)
        skip("py_outrider not available for testing")
}

test_that("Test findEncodingDim PROTRIDER", {
    
    # Test Find Encoding Dimension.
    skip_if_no_py_outrider()
    ods <- makeExampleProtriderDataSet(40,40)
    metadata(ods)$encDimTable <- NULL
    metadata(ods)$optimalEncDim <- NULL
    ods <- filterExpression(ods)
    intensbefore <- observed(ods)
    expect_equal(getBestQ(ods), NA_integer_)
    ods <- findEncodingDim(ods, params=3:5, iterations=2, usePython=TRUE)
    
    expect_equal(intensbefore, observed(ods))
    expect_is(ods, 'Outrider2DataSet')
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
    tCts <- assay(ods, 'trueObservations') # tCts <- assay(ods, 'trueCounts')
    numExpect <- nrow(ods) * ncol(ods) * freq
    
    expect_equal(counts(ods)[tCor == 0], tCts[tCor == 0])
    expect_true(all(counts(ods)[tCor != 0] != tCts[tCor != 0]))
})
