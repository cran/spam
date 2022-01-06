# HEADER ####################################################
# This is file spam/tests/testthat/test-definitions.R.      #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################
rm(list = ls())
source("helper.R")

## library("testthat")
## library("spam64", lib.loc = LIB.LOC)
## library("spam", lib.loc = "../../../lib/")


context("test-definitions.R")


################## PR 18272 #################
# https://bugs.r-project.org/show_bug.cgi?id=18272
# Patch by Mike Chirico

Sys.setenv("_R_CHECK_LENGTH_1_LOGIC2_" = "true")
# all.equal(c(1, 1), c(1.01, 1.01), scale = c(.01, .01))  # patch 2022.01.03
all.equal.spam(diag.spam(2), 1.01*diag.spam(2), scale = c(.01, .01))  # patch 2022.01.03


Sys.unsetenv("_R_CHECK_LENGTH_1_LOGIC2_")
# all.equal(c(1, 1), c(1.01, 1.01), scale = c(.01, .01)) # patch 2022.01.03
all.equal.spam(diag.spam(2), 1.01*diag.spam(2), scale = c(.01, .01))  # patch 2022.01.03




###################################
# general stuff, mainly from help of all.equal.spam
test_that("all.equal.spam", {
  obj <- diag.spam(2)
  obj[1,2] <- .Machine$double.eps

  expect_equivalent( all.equal( diag.spam(2), obj),"Lengths (2, 3) differ")
  expect_equivalent( all.equal( t(obj), obj), c("Column-sparsity structure differ (at least 1 instance(s))",
      "Row-sparsity structure differ (at least 1 instance(s))"))
  expect_equivalent( all.equal( t(obj), obj*1.1), c("Column-sparsity structure differ (at least 1 instance(s))",
      "Row-sparsity structure differ (at least 1 instance(s))", "Mean relative difference: 0.1" ))

  # We can compare a spam to a matrix
  expect_equivalent( all.equal(diag(2),diag.spam(2)), TRUE)


 })

