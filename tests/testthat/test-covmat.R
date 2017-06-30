# This is file ../spam/tests/covmat.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     
rm(list = ls())
source("helper.R")

## library("testthat")
## library("spam64", lib.loc = LIB.LOC)
## library("spam", lib.loc = "../../../lib/")


context("test-covmat.R")



test_that("cov.*", {
    h <- nearest.dist(100*1:10, 100*1:10+1:10, delta=10)
    expect_identical(cov.exp(1:10, 10), cov.exp(h, 10)@entries)
    expect_identical(cov.sph(1:10, 10), cov.sph(h, 10)@entries)
    expect_identical(cov.nug(1:10, 10), cov.nug(h, 10)@entries)
    expect_identical(cov.wu1(1:10, 10), cov.wu1(h, 10)@entries)
    expect_identical(cov.wu2(1:10, 10), cov.wu2(h, 10)@entries)
    expect_identical(cov.wu3(1:10, 10), cov.wu3(h, 10)@entries)
    expect_identical(cov.wend1(1:10, 10), cov.wend1(h, 10)@entries)
    expect_identical(cov.wend2(1:10, 10), cov.wend2(h, 10)@entries)
    expect_identical(cov.mat(1:10, 10), cov.mat(h, 10)@entries)
})
