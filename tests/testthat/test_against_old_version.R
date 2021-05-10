context("Testing against old version: ")

test_that("test against old version", {
    odsOld <- get_test_outrider_dataset()
    
    ods <- OutriderDataSet(countData=counts(odsOld))
    sizeFactors(ods) <- sizeFactors(odsOld)
    normalizationFactors(ods) <- normalizationFactors(odsOld)
    
    ods <- fit(ods)
    expect_equal(round(mcols(ods)[['mu']], 4), mcols(odsOld)[['mu']])
    expect_equal(round(theta(ods), 4), mcols(odsOld)[['disp']])
    
    ods <- preprocess(ods)
    ods <- computeZscores(ods)
    #expect_equal(round(zScore(ods), 1), round(zScore(odsOld), 1))
    
    ods <- computePvalues(ods, method='BH')
    expect_equal(round(pValue(ods), 4), assay(odsOld, 'pValue'))
    expect_equal(round(padj(ods), 4), assay(odsOld, 'padjust'))
})
