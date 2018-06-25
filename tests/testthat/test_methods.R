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
    ods <- makeExampleOutriderDataSet()
    expect_error(results(ods), "The P-values are not computed yet..*")
})

test_that("normalization method", {
    ods <- makeExampleOutriderDataSet()[seq_len(3),1]
    counts(ods) <- matrix(c(1L,2L,3L), ncol=1)
    expect_null(normalizationFactors(ods))
    nMat <- matrix(5L, ncol=1, nrow=3)
    normalizationFactors(ods) <- nMat
    expect_equivalent(normalizationFactors(ods), nMat)
    expect_error(results(ods), "The P-values are not computed yet..*")
})

test_that("fit method", {
    ods <- makeExampleOutriderDataSet(50)
    ods <- estimateSizeFactors(ods)
    ods <- fit(ods)
    expect_is(mcols(ods)[,"mu"], "numeric")
    expect_equal(length(mcols(ods)[,"mu"]), nrow(ods))
    expect_is(mcols(ods)[,"disp"], "numeric")
    expect_equal(length(mcols(ods)[,"disp"]), nrow(ods))
})
