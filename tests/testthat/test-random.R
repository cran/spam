# HEADER ####################################################
# This is file spam/tests/testthat/test-random.R.           #
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

context("test-random.R")

test_that("check arguments", {
  rspam1 <- spam::spam_random(7, digits = 2, distribution = runif, min = -1, sym = TRUE)
  expect_identical(spam::isSymmetric.spam(rspam1), TRUE)

  rspam2 <- spam::spam_random(7, digits = 2, distribution = NULL, sym = FALSE, spd = TRUE)
  cholrspam2 <- chol(rspam1)

  expect_identical(spam::isSymmetric.spam(rspam2), TRUE)
  expect_equivalent(class(cholrspam2), "spam.chol.NgPeyton")
})
