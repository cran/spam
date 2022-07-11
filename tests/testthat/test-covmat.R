# HEADER ####################################################
# This is file spam/tests/testthat/test-covmat.R.           #
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


context("test-covmat.R")


test_that("cov.sph", {
    nl <- 10
    range <- .5
    h <- nearest.dist(cbind( runif(nl), runif(nl)),  delta=range, upper=NULL)
    expect_identical( cov.sph(h, c(range, 1, 0)),   cor.sph(h, range) )

    powerboost()
    expect_identical( cov.sph(h, c(range, 1+1e-10 )),     cov.sph(h, range) )
    expect_identical( cov.sph(h, c(range, 1+1e-10 )),     cov.sph(h, range) )
    expect_identical( cov.sph(h, c(range, 1+1e-9, 1e-9)), cov.sph(h, range) )

    powerboost("aus")    
    expect_equal( cov.sph(h, c(range, 1+1e-9, 1e-9)), cov.sph(h, range) )
})


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

test_that("zero range", {
  set.seed(42)
  n <- 100
  sampleData <- matrix(rnorm(n), nrow = sqrt(n))
  sampleData[4:6, 6:8] <- sampleData[4:6, 6:8] + 5

  x <- 1:10
  locs <- expand.grid(x, x)
  dim1 <- dim2 <- 10
  distmat <- nearest.dist( locs, upper=NULL, delta = 1000) # distance matrix
  out <- spam::mle.nomean.spam(y = (c(sampleData)), distmat = distmat, Covariance =
                          cov.mat,
                        theta0 = c(2, 1, 1, .01),
                        thetalower = c(0, 0.1, 0.1, 0.1),
                        thetaupper = c(5, 5, 5, .2))
})




test_that("cov.matXY*", {
    h <- seq(0, to=10, len=150)
    theta <- c(1.4, 1, 1)
    theta1 <- c(1.4, 1, 0.5, 1)
    theta3 <- c(1.4, 1, 1.5, 1)
    theta5 <- c(1.4, 1, 2.5, 1)
    
    expect_equal(cov.exp(h, theta), cov.mat(h,theta1))
    expect_equal(cov.mat12(h, theta), cov.mat(h,theta1))
    expect_equal(cov.mat32(h, theta), cov.mat(h,theta3))
    expect_equal(cov.mat52(h, theta), cov.mat(h,theta5))
})


theta <- c(1.4, 1, 1)
theta1 <- c(1.4, 1, 0.5, 1)
theta3 <- c(1.4, 1, 1.5, 1)
theta5 <- c(1.4, 1, 2.5, 1)
# h  <- seq(0, to=10, len=10000)
# microbenchmark::microbenchmark(a=cov.exp(h, theta), a2=cov.mat(h,theta1), b=cov.mat32(h, theta), b2=cov.mat(h,theta3), c=cov.mat52(h, theta), c2=cov.mat(h,theta5))
