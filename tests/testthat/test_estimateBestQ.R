context("Testing the estimateBestQ function")

library(denoiseR)

# Hyperparameter optimization (grid-search)

test_that("Test hyperparameter optimization", {
  
  ods <- makeExampleOutriderDataSet(dataset = 'Kremer')
  ods <- ods[1:40, 1:40]
  ods <- filterExpression(ods, filterGenes=TRUE, minCounts=TRUE)
  countsbefore <- counts(ods)
  expect_equal(getBestQ(ods), NA_integer_)
  ods <- estimateBestQ(ods, useOHT=FALSE, params=3:5, iteration=2)
  
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

# Optimal Hard Thresholding

test_that("Input validation handles NULL and non-matrix inputs", {
  expect_error(estimateBestQ(), 
               "Please provide an OutriderDataSet or a z-score matrix.")
  expect_error(estimateBestQ(NULL, NULL), 
               "Please provide an OutriderDataSet or a z-score matrix.")
  expect_error(estimateBestQ(zScoresOHT = "not a matrix"), 
               "Provided zScoresOHT are not a matrix.")
  expect_error(estimateBestQ(ods = "not an ods"), 
               "Please provide an OutriderDataSet.")
  
  ctsFile <- system.file('extdata', 'GTExSkinSmall.tsv',
                         package='OUTRIDER')
  ctsTable <- read.table(ctsFile, check.names=FALSE)
  ods <- OutriderDataSet(countData=ctsTable)
  # filter out non expressed genes
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  ods <- estimateSizeFactors(ods)
  
  expect_warning(estimateBestQ(ods = ods, zScoresOHT = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)),
                 "Provided z-scores are ignored and recalculated from ods.")
})

test_that("User is notified about invalid matrix values", {
  expect_error(estimateBestQ(zScoresOHT = matrix(c(1, 2, 3, 4, 5, Inf), nrow = 3, ncol = 2)), 
               "Z-score matrix contains infinite values.")
})

test_that("optimalSVHTCoef works correctly", {
  expect_equal(optimalSVHTCoef(0.5),
               1.98, tolerance = 0.01)
  expect_equal(optimalSVHTCoef(0.1),
               1.58, tolerance = 0.01)
})

test_that("medianMarchenkoPastur works correctly", {
  # Expected outputs are derived from Table IV in Gavish and Donoho (2014)
  expect_equal(optimalSVHTCoef(0.5) / sqrt(medianMarchenkoPastur(100, 200)),
               2.1711, tolerance = 0.0001)
  expect_equal(optimalSVHTCoef(0.1) / sqrt(medianMarchenkoPastur(100, 1000)),
               1.6089, tolerance = 0.0001)
})

test_that("Encoding dimensions are properly calculated for simulated z-scores", {
  # Simulate zScore matrix consisting of signal and noise
  set.seed(42)
  numGenes <- 10000
  numSamples <- 200
  latentDim <- 50
  signalNoiseRatio <- 5
  zTilde <- LRsim(numGenes, numSamples, latentDim, signalNoiseRatio)$X * 1000
  
  expect_equal(estimateBestQ(zScoresOHT = zTilde),
               latentDim)
  
  # Simulate zScore matrix with beta > 1
  set.seed(42)
  numGenes <- 50
  numSamples <- 200
  latentDim <- 20
  signalNoiseRatio <- 5
  zTilde <- LRsim(numGenes, numSamples, latentDim, signalNoiseRatio)$X * 1000
  
  expect_equal(estimateBestQ(zScoresOHT = zTilde),
               latentDim)
  
  # Simulate zScore matrix consisting of noise only
  set.seed(42)
  latentDim <- 0
  zTilde <- matrix(rnorm(numGenes * numSamples), nrow = numGenes, ncol = numSamples)
  expect_warning(expect_equal(estimateBestQ(zScoresOHT = zTilde), 2),
               paste("Optimal latent space dimension is smaller than 2\\. Check your count matrix and",
                 "verify that all samples have the expected number of counts",
                 "\\(hist\\(colSums\\(counts\\(ods\\)\\)\\)\\)\\.",
                 "For now\\, the latent space dimension is set to 2\\.", collapse = "\n"))
})

test_that("Encoding dimensions are properly calculated for real ODS", {
  ctsFile <- system.file('extdata', 'GTExSkinSmall.tsv',
                         package='OUTRIDER')
  ctsTable <- read.table(ctsFile, check.names=FALSE)
  ods <- OutriderDataSet(countData=ctsTable)
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  ods <- estimateSizeFactors(ods)
  
  outsingleResult <- 5 # Expected value was calculated with OutSingle
  ods <- estimateBestQ(ods=ods, useOHT=TRUE)
  
  expect_equal(metadata(ods)[['optimalEncDim']], getBestQ(ods))
  expect_equal(metadata(ods)[['optimalEncDim']], outsingleResult,
               tolerance = 1)
})