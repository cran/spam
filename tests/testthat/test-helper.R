# HEADER ####################################################
# This is file spam/tests/testthat/test-helper.R.           #
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


context("test-helper.R")



########


set.seed(14)

# bdiag.spam:

A <- spam(rnorm(10),2)
B <- spam(rnorm(16),4)

test_that("bdiag.spam", {
    spamtest_eq( bdiag.spam(A),A)

    spamtest_eq( bdiag.spam(A,B),rbind(cbind(A,rep(0,8)),
                                       cbind(spam(rep(0,20),4),B)))
})


