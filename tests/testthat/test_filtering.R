context("Testing filtering: ")

test_that("Expression filtering", {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    map <- select(org.Hs.eg.db, keys=keys(txdb, keytype = "GENEID"), 
                  keytype="ENTREZID", columns=c("SYMBOL"))
    
    ods <- makeExampleOutriderDataSet(dataset='KremerNBader')
    
    expect_error(filterExpression(ods), 
            'rowRanges\\(object\\) has all ranges of.* FPKM values')
    
    expect_warning(filterExpression(ods, gtfFile=txdb), 
            'Some genes .*999.* are not found')
    
    expect_warning(filterExpression(ods, gtfFile=txdb, mapping=map, save=TRUE), 
                   'Some genes .*n=17[72]\\).* are not found')
    
})
