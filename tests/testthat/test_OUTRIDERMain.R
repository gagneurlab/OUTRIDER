context("Testing OUTRIDER main function")

test_that("test main OUTRIDER function", {
    maxIter <- 6
    ods <- makeExampleOutriderDataSet(30, 30)
    ods <- preprocess(ods) # ods <- estimateSizeFactors(ods)
    ods <- controlForConfounders(ods, iteration=maxIter)
    ods <- computePvalues(ods)
    # ods <- computeZscores(ods)
    ods <- computeEffectSizes(ods, distribution="nb", 
                                effect_types=c("fold_change", "zscores"))
    res <- results(ods, all=TRUE)
    expect_is(res, "data.table")
    
    ods2 <- OUTRIDER(ods, iteration=maxIter, implementation="autoenc")
    ods3 <- OUTRIDER(ods, implementation="pca")
    res2 <- results(ods2, all=TRUE)
    
    expect_equal(ods, ods2)
    expect_equal(res, res2)
    expect_is(ods3, "OutriderDataSet")
    expect_is(pValue(ods3), "matrix")
    expect_equal(dim(pValue(ods3)), dim(pValue(ods2)))
})

# helper function to skip tests if we don't have the 'py_outrider' module
skip_if_no_py_outrider <- function() {
    have_py_outrider <- reticulate::py_module_available("py_outrider")
    if (!have_py_outrider)
        skip("py_outrider not available for testing")
}

test_that("test main OUTRIDER function for PROTRIDER", {
    skip_if_no_py_outrider()
    maxIter <- 6
    ods <- makeExampleProtriderDataSet(30, 30)
    ods <- preprocess(ods) # ods <- estimateSizeFactors(ods)
    ods <- controlForConfounders(ods, iterations=maxIter)
    ods <- computePvalues(ods, distribution="gaussian")
    ods <- computeEffectSizes(ods, distribution="gaussian", 
                              effect_types=c("fold_change", "zscores"))
    res <- results(ods, all=TRUE)
    expect_is(res, "data.table")
    
    ods2 <- OUTRIDER(ods, iterations=maxIter, latent_space_model="autoenc", 
                        decoder_model="autoenc")
    ods3 <- OUTRIDER(ods, latent_space_model="pca", decoder_model="pca", 
                        iterations=1)
    res2 <- results(ods2, all=TRUE)
    
    expect_equal(ods, ods2)
    expect_equal(res, res2)
    expect_is(ods3, "Outrider2DataSet")
    expect_is(pValue(ods3), "matrix")
    expect_equal(dim(pValue(ods3)), dim(pValue(ods2)))
})
