Sys.setenv("R_TESTS"="")
library(testthat)
library(OUTRIDER)

register(SerialParam())

test_check("OUTRIDER")
