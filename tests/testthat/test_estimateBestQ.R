context("Testing the estimateBestQ function (Optimal Hard Thresholding)")

test_that("Input validation handles NULL and non-matrix inputs", {
  expect_error(findEncodingDimOptimalHardThreshold(), 
               "Please provide an OutriderDataSet or a z-score matrix.")
  expect_error(findEncodingDimOptimalHardThreshold(NULL, NULL), 
               "Please provide an OutriderDataSet or a z-score matrix.")
  expect_error(findEncodingDimOptimalHardThreshold(zScores = "not a matrix"), 
               "Provided zScores are not a matrix.")
  expect_error(findEncodingDimOptimalHardThreshold(ods = "not an ods"), 
               "Please provide an OutriderDataSet.")
  
  ctsFile <- system.file('extdata', 'GTExSkinSmall.tsv',
                         package='OUTRIDER')
  ctsTable <- read.table(ctsFile, check.names=FALSE)
  ods <- OutriderDataSet(countData=ctsTable)
  # filter out non expressed genes
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  
  expect_warning(findEncodingDimOptimalHardThreshold(ods = ods, zScores = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)),
                 "Provided z-scores are ignored and recalculated from ods.")
})

test_that("User is notified about invalid matrix dimensions or values", {
  expect_error(findEncodingDimOptimalHardThreshold(zScores = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), 
                                                                    nrow = 2, ncol = 4)),
               c("Number of columns (samples) is larger than number of rows (genes).",
                 "OHT does not work for such cases."))
  expect_error(findEncodingDimOptimalHardThreshold(zScores = matrix(c(1, 2, 3, 4, 5, Inf), nrow = 3, ncol = 2)), 
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
  
  expect_equal(findEncodingDimOptimalHardThreshold(zScores = zTilde),
               latentDim)
  
  # Simulate zScore matrix consisting of noise only
  set.seed(42)
  latentDim <- 0
  zTilde <- matrix(rnorm(numGenes * numSamples), nrow = numGenes, ncol = numSamples)
  expect_error(expect_equal(findEncodingDimOptimalHardThreshold(zScores = zTilde), 
                            latentDim),
               c("Latent space dimension is smaller than 2. Check your count matrix and",
                 "verify that all samples have the expected number of counts.",
                 "hist(colSums(counts(ods)))"))
})

test_that("Encoding dimensions are properly calculated for real ODS", {
  ctsFile <- system.file('extdata', 'GTExSkinSmall.tsv',
                         package='OUTRIDER')
  ctsTable <- read.table(ctsFile, check.names=FALSE)
  ods <- OutriderDataSet(countData=ctsTable)
  ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
  outsingleResult <- 5 # Expected value was calculated with OutSingle
  expect_equal(findEncodingDimOptimalHardThreshold(ods = ods), outsingleResult, 
               tolerance = 1)
})