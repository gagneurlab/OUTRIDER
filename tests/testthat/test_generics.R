context("Testing generic functions")

test_that("setter and getter function", {
    ord <- makeExampleOutriderDataSet(30, 30)
    sf <- runif(dim(ord)[2], 0.8, 1.2)
    sizeFactors(ord) <- sf
    
    expect_equal(sizeFactors(ord), sf, check.names=FALSE)
    expect_equal(t(t(counts(ord))/sf), counts(ord, normalized=TRUE))
})

test_that("sizefactor estimation", {
    cts <- getCountDataWithNames(matrix(rep(c(1:10),4), nrow=10) + 1:40)
    ord <- OutriderDataSet(countData=cts)
    exp_sf <- c(0.471, 0.905, 1.333, 1.763)
    sf <- sizeFactors(estimateSizeFactors(ord))
    
    expect_equal(sf, exp_sf, tolerance=.001, check.names=FALSE)
})

test_that("normalization function", {
    n <- 30
    m <- 4
    ord <- makeExampleOutriderDataSet(n=n, m=m)
    
    normF <- matrix(c(1:m)/m, nrow=n, ncol=m)
    normalizationFactors(ord) <- normF
    
    expect_equivalent(normF, normalizationFactors(ord))
    E <- t(apply(normF, 1, pmax, 0.5))
    expect_equal(counts(ord)/E * exp(rowMeans(log(E))), 
            counts(ord, normalized=TRUE))
    
    # check if subsetting affects the computation of the log geom mean
    expect_equal((counts(ord)/E * exp(rowMeans(log(E))))[1:10,1:3], 
            counts(ord[1:10,1:3], normalized=TRUE))
    
    normalizationFactors(ord) <- NULL
    expect_null(normalizationFactors(ord))
    
    ord <- OutriderDataSet(countData=getCountDataWithNames(
            matrix(c(1,1,3,3), ncol=2)))
    colData(ord)[['sizeFactor']] <- c(.5, 1.5)
    expect_true(all(
            counts(ord, normalized = TRUE) == matrix(c(2,2,2,2), ncol = 2)))
    
    ord <- OutriderDataSet(countData=getCountDataWithNames(
            matrix(c(1,1,3,3), ncol=2)))
    normalizationFactors(ord) <- matrix(c(.5, .5, 1.5, 1.5), ncol=2)
    expect_true(all(round(counts(ord, normalized=TRUE), 4) == 
            matrix(1.7321, ncol=2, nrow=2)))
})

test_that("counts_function", {
    ord <- makeExampleOutriderDataSet(30, 30)
    cts <- counts(ord)
    counts(ord) <- cts+1L
    
    expect_equal(assays(ord)[['counts']], cts+1L)
    expect_equal(counts(ord), assays(ord)[['counts']])
})

