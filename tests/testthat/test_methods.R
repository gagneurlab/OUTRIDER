context("Testing computational functions")

test_that("filter expression", {
    ods <- makeExampleOutriderDataSet(dataset="GTExSkinSmall")
    annotationFile <- system.file("extdata", "gencode.v19.genes.small.gtf.gz",
                                  package="OUTRIDER")
    ods <- filterExpression(ods, annotationFile, filterGenes=TRUE,
                               fpkmCutoff=0)
    
    expect_true(all(mcols(ods)[,"passedFilter"]))
    expect_equal(dim(ods), c(82, 50))
    expect_true("basepairs" %in% colnames(mcols(ods)))
})

test_that("fitting method", {
    ods <- makeExampleOutriderDataSet()
    expect_error(fit(ods), "Please provide sizeFactors or normal.*")
})

test_that("pvalue calculation", {
    ods <- makeExampleOutriderDataSet()
    expect_error(computePvalues(ods), "Please fit the models first to .*")
})

test_that("result method", {
    ods <- makeExampleOutriderDataSet(150, 50)
    expect_error(results(ods), "Please calculate..*")
    ods <- OUTRIDER(ods, q=2, iteration=2)
    
    expect_warning(res <- results(ods, padj=1e-15), "No significant events:")
    expect_equal(colnames(res), colnames(results(ods)))
    res_all <- results(ods, all=TRUE)
    expect_equal(nrow(res_all), nrow(ods)*ncol(ods))
    expect_true(all(results(ods)$aberrant))
    expect_equal(nrow(results(ods, padjCutoff=0.05, zScoreCutoff=0, round=TRUE)), 
        nrow(res_all[padjust <= 0.05 & abs(zScore) >= 0])
    )
})

test_that("normalization method", {
    ods <- makeExampleOutriderDataSet(5, 5)
    counts(ods) <- matrix(1:5, ncol=5, nrow=5, dimnames=dimnames(ods))
    expect_null(normalizationFactors(ods))
    nMat <- matrix(6:10, ncol=5, nrow=5)
    normalizationFactors(ods) <- nMat
    expect_equivalent(normalizationFactors(ods), nMat)
    expect_error(results(ods), "Please calculate..*")
})

test_that("fit method", {
    ods <- makeExampleOutriderDataSet(30, 30)
    ods <- estimateSizeFactors(ods)
    ods <- fit(ods)
    expect_is(mcols(ods)[,"mu"], "numeric")
    expect_equal(length(mcols(ods)[,"mu"]), nrow(ods))
    expect_is(theta(ods), "numeric")
    expect_equal(length(theta(ods)), nrow(ods))
})
