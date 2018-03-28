context("Testing computational functions: ")

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
    ods <- makeExampleOutriderDataSet()
    counts(ods[1:3,1]) <- matrix(c(1L,2L,3L), ncol = 1)
    normalizationFactors(ods[1:3,1]) <- matrix(c(5L,5L,5L), ncol = 1)
    expect_error(results(ods), "The P-values are not computed yet..*")
    
})

test_that("fit method", {
    ods <- makeExampleOutriderDataSet(50)
    ods <- estimateSizeFactors(ods)
    ods <- fit(ods)
    expect_true(all(mcols(ods)[,"mu"]))
    expect_true(all(mcols(ods)[,"disp"]))
})

test_that("test math", {
    # raw data
    cts <- matrix(byrow=TRUE, nrow=3, as.integer(c(
            c( 4,  1,  8,  1,  9,  3,  0,  1,  1,  3,  1,  9, 25,  6,  1),
            c(15, 21, 18,  8, 16, 25, 12, 22, 28,  6, 10, 31, 10, 11,  6),
            c(15, 17, 32, 10, 23, 32,  8,  6, 58, 20, 11, 12,  8, 24, 12))))
    sizeF <- c(0.925, 0.974, 0.889, 0.950, 0.943, 1.010, 1.023, 1.130, 1.100, 
            0.999, 0.866, 1.152, 0.805, 1.174, 1.120)
    
    # expected data
    mu_round4 <- c(5.1835, 15.7957, 19.0995)
    disp_round4 <- c(0.8608, 6.6356, 3.3016)
    pval_round4 <- matrix(byrow=TRUE, nrow=3, c(
            c(0.9796, 0.9916, 0.5348, 0.6722, 0.8275, 1.0000, 0.5510, 0.5971, 
                    0.6083, 1.0000, 0.7147, 0.5024, 0.0399, 0.7912, 0.6008),
            c(0.9796, 0.9916, 0.5348, 0.6722, 0.8275, 0.4178, 0.6948, 0.5971,
                    0.3242, 0.4463, 0.7147, 0.4869, 0.8042, 0.7912, 0.3026), 
            c(0.9796, 1.0000, 0.5348, 0.6722, 0.8275, 0.4178, 0.5510, 0.4554,
                    0.0744, 1.0000, 0.7147, 0.5024, 0.7466, 0.7912, 0.6008)))
    zscore_round3 <- matrix(byrow=TRUE, nrow=3, c(
            c( 0.340, -0.701,  1.014, -0.675,  1.064,  0.005, -1.500, -0.861, 
                    -0.832,  0.017, -0.575,  0.849,  2.262,  0.445, -0.852),
            c( 0.269,  0.828,  0.712, -0.991,  0.356,  1.101, -0.376,  0.610,
                    +1.151, -1.623, -0.377,  1.261, -0.225, -0.832, -1.862),
            c( 0.023,  0.137,  1.334, -0.666,  0.686,  1.115, -1.137, -1.739),
                    +1.965,  0.358, -0.357, -0.710, -0.726,  0.380, -0.661))
    
    # run deterministic pipeline
    ods <- OutriderDataSet(countData=cts)
    sizeFactors(ods) <- sizeF
    ods <- fit(ods)
    ods <- computePvalues(ods)
    ods <- computeZscores(ods)
    
    # test results
    expect_equal(round(mcols(ods)[['mu']], 4), mu_round4)
    expect_equal(round(mcols(ods)[['disp']], 4), disp_round4)
    expect_equivalent(round(assays(ods)[['padjust']], 4), pval_round4)
    expect_equivalent(round(assays(ods)[['zScore']], 3), zscore_round3)
})