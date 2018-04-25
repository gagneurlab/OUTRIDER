context("Testing helper functions: ")

test_that("read/write NB model", {
    ods <- makeExampleOutriderDataSet(n=30, m=30)
    rown <- rownames(ods)
    filename <- "test_nb_model.txt"
    
    ###
    ### writing model
    ###
    expect_error(writeNBModel(ods, NA), 'Please provide a correct file name .*')
    expect_error(writeNBModel(ods, 1L), 'Please provide a correct file name .*')
    
    rownames(ods) <- NULL
    expect_error(writeNBModel(ods, filename),
            'Please provide rownames, if you want to save .*')
    
    rownames(ods) <- rown[c(1:15,1:15)]
    expect_error(writeNBModel(ods, filename),
                 'Please provide uniq rownames, if you want to save .*')
    
    rownames(ods) <- rown
    expect_error(writeNBModel(ods, filename),
                 'The following mcols are missing: .* Please run first.*')
    
    ods <- estimateSizeFactors(ods)
    ods <- fit(ods)
    expect_null(writeNBModel(ods, filename))
    expect_true(file.exists(filename))
    
    
    ###
    ### reading model
    ###
    expect_error(readNBModel(ods, NULL), "Please provide a correct file name.*")
    expect_error(readNBModel(ods, 1L), "Please provide a correct file name.*")
    expect_error(readNBModel(ods, NA), "Please provide a correct file name.*")
    expect_error(readNBModel(ods, "NULL"), "Can not find given model.*")
    
    rownames(ods) <- NULL
    expect_error(readNBModel(ods, filename), 'Please provide rownames to.*')
    
    ods <- makeExampleOutriderDataSet(n=30, m=30)
    ods <- readNBModel(ods, filename)
    expect_true(all(c('loggeomeans', 'mu', 'disp') %in% colnames(mcols(ods))))
})
