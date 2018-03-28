context("Testing against old version: ")

test_that("test against old version", {
    odsOld <- get_test_outrider_dataset()
    
    ods <- OutriderDataSet(countData=counts(odsOld))
    sizeFactors(ods) <- sizeFactors(odsOld)
    normalizationFactors(ods) <- normalizationFactors(odsOld)
    
    ods <- fit(ods)
    expect_equal(round(mcols(ods)[['mu']], 4), mcols(odsOld)[['mu']])
    expect_equal(round(mcols(ods)[['disp']], 4), mcols(odsOld)[['disp']])
    
    ods <- computeZscores(ods)
    expect_equal(round(assay(ods, 'zScore'), 4), assays(odsOld, 'zScore'))
    
    ods <- computePvalues(ods)
    expect_equal(round(assay(ods, 'pValue'), 4), assay(odsOld, 'pValue'))
    expect_equal(round(assay(ods, 'padjust'), 4), assay(odsOld, 'padjust'))
})