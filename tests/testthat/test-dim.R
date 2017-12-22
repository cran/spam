# HEADER ####################################################
# This is file spam/tests/testthat/test-dim.R.              #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      #
#    Daniel Gerber [ctb], Kaspar Moesinger [ctb],           #
#    Youcef Saad [ctb] (SPARSEKIT),                         #
#    Esmond G. Ng [ctb] (Fortran Cholesky routines),        #
#    Barry W. Peyton [ctb] (Fortran Cholesky routines),     #
#    Joseph W.H. Liu [ctb] (Fortran Cholesky routines),     #
#    Alan D. George [ctb] (Fortran Cholesky routines),      #
#    Esmond G. Ng [ctb] (Fortran Cholesky routines),        #
#    Barry W. Peyton [ctb] (Fortran Cholesky routines),     #
#    Joseph W.H. Liu [ctb] (Fortran Cholesky routines),     #
#    Alan D. George [ctb] (Fortran Cholesky routines)       #
# HEADER END ################################################

rm(list = ls())
source("helper.R")

## library("testthat")
## library("spam64", lib.loc = LIB.LOC)
## library("spam", lib.loc = "../../../lib/")


context("test-dim.R")


# simple tests:
########################################################################


# construct matrices:
n <- 10
m <- 15

set.seed(14)
tt <- matrix(rnorm(m*n),n,m)
tt[tt<0] <- 0

ss <- as.spam(tt)



test_that("dim", {
    spamtest_eq(ss,tt)
    
    dim(ss) <- c(m,n)
    dim(tt) <- c(m,n)
    spamtest_eq(ss,tt)
    
    dim(ss) <- c(m*n,1)
    dim(tt) <- c(m*n,1)
    spamtest_eq(ss,tt)
    
    dim(ss) <- c(1, m*n)
    dim(tt) <- c(1, m*n)
    spamtest_eq(ss,tt)
    
    expect_error( dim(ss) <- c(-1, -m*n),
                 "Indices need to be positive")
    expect_error( dim(ss) <- c(1, m, n),
                 "dims should be of length 1 or 2")
})
