

# Test Find Encoding Dimension.
ods <- makeExampleOutriderDataSet(dataset = 'Kremer')
res <- findEncodingDim(ods)
is(res, 'data.frame')

# Test In silico outliers.
freq <- 1E-2
zScore=3
inj='both'
ods <- injectOutliers(ods, freq=freq, zScore=zScore, inj=inj)
#Check that everything worked
expect_true(all(counts(ods)[
    assays(ods)[['trueCorruptions']]==0] == assay(ods, 'trueCounts')[
        assays(ods)[['trueCorruptions']]==0]))
expect_true(any(counts(ods)[
    assays(ods)[['trueCorruptions']]!=0] != assay(ods, 'trueCounts')[
        assays(ods)[['trueCorruptions']]!=0]))

