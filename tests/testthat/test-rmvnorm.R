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


context("test-rmvnorm.R")



########
verbose <- FALSE

set.seed(14)


n <- 5
Sigma <- .25^abs(outer(1:n,1:n,'-'))
Q <- as.spam(solve(Sigma))
b <- 1:n

struct <- chol(Q, verbose=verbose)

test_that("rmvnorm", {

    set.seed(14)
    tmp1 <- rmvnorm(10,  Sigma=Sigma, verbose=verbose)
    set.seed(14)
    tmp2 <- rmvnorm.spam(10, Sigma=Sigma)
    expect_equal( tmp1, tmp2)

    set.seed(14)
    tmp1 <- rmvnorm(10,  Sigma=Q)
    set.seed(14)
    tmp2 <- rmvnorm.spam(10, Sigma=Q)
    expect_equal( tmp1, tmp2)



    set.seed(14)
    tmp1 <- rmvnorm.canonical(10, b, Q)
    set.seed(14)
    spamtest_eq( rmvnorm.canonical(10, b, Q, Lstruct=struct), tmp1 )
    set.seed(14)
    spamtest_eq( rmvnorm.prec(10, solve(Q,b), Q), tmp1 )
    set.seed(14)
    spamtest_eq( rmvnorm.prec(10, solve(Q,b), Q, Lstruct=struct), tmp1 )
    set.seed(14)
    ## cat("For rmvnorm.canonical:\n- comparing sample mean with truth:\n")
    ## for (i in 10^(1:4))
    ##     cat('    sample size n=',i,' yields  Frobenius-norm:',
    ##         norm( apply(rmvnorm.canonical(i, b, Q, Lstruct=struct), 2,mean)- solve(Q,b),'f'),'\n')
    ## cat("- comparing sample variance with truth:\n")
    ## for (i in 10^(1:4)){
    ##     cat('    sample size n=',i,' yields Frobenius-norm:',
    ##         norm( var( rmvnorm.canonical(i, b, Q=Q, Lstruct=struct))- Sigma,'f'),'\n')
    ## set.seed(14)
    ## cat("For rmvnorm.prec:\n- comparing sample mean with truth:\n")
    ## for (i in 10^(1:4))
    ##     cat('    sample size n=',i,' yields  Frobenius-norm:',
    ##         norm( apply(rmvnorm.prec(i, b, Q, Lstruct=struct), 2,mean)- b,'f'),'\n')
    ## cat("- comparing sample variance with truth:\n")
    ## for (i in 10^(1:4)){
    ##     cat('    sample size n=',i,' yields Frobenius-norm:',
    ##         norm( var( rmvnorm.prec(i, Q=Q, Lstruct=struct))- Sigma,'f'),'\n')
    ## }
    ## set.seed(14)
    ## cat("For rmvnorm.spam:\n- comparing sample mean with truth:\n")
    ## for (i in 10^(1:4))
    ##   cat('    sample size n=',i,' yields  Frobenius-norm:',
    ##       norm( apply(rmvnorm.spam(i, b, as.spam(Sigma), Lstruct=struct), 2,mean)- b,'f'),'\n')
    ## cat("- comparing sample variance with truth:\n")
    ## for (i in 10^(1:4)){
    ##   cat('    sample size n=',i,' yields Frobenius-norm:',
    ##       norm( var( rmvnorm.spam(i, b, as.spam(Sigma), Lstruct=struct))- Sigma,'f'),'\n')
    ## }
})


test_that("rmvt", {

  set.seed(14)
  tmp1 <- rmvt.spam(10,  Sigma=Sigma, df = 0)
  set.seed(14)
  tmp2 <- rmvnorm.spam(10, Sigma = Sigma)
  expect_equal( tmp1, tmp2)

  set.seed(14)
  tmp1 <- rmvt.spam(10, Sigma = Sigma, df = 4)
  set.seed(14)
  tmp2 <- rmvnorm.spam(10, Sigma = Sigma)/sqrt(rchisq(10, 4)/4)
  expect_equal( tmp1, tmp2)

  set.seed(14)
  tmp1 <- rmvt.spam(10, Sigma = Sigma, df = 4, type = "shifted")
  set.seed(14)
  tmp2 <- sweep(rmvnorm.spam(10, Sigma = Sigma)/sqrt(rchisq(10, df = 4)/4), 2, rep(0, nrow(Sigma)), "+")
  expect_equal( tmp1, tmp2)

  set.seed(14)
  tmp1 <- rmvt(10, Sigma = Sigma, df = 4, type = "shifted")
  set.seed(14)
  tmp2 <- rmvt.spam(10, Sigma = Sigma, df = 4, type = "shifted")
  expect_equal( tmp1, tmp2)
})



require('fields')

context("test-rgrf.R: locations")


tau <- 0
nx <- 49
ny <- nx
field <- rgrf(1, nx=nx, ny=ny, tau=tau, Covariance="cov.nug", theta=1)
quilt.plot(cbind(attr(field,"locs"),z=field), nx=nx, ny=ny)
points(attr(field,"locs"))

expect_error(rgrf(1, nx=3, tau=.5, Covariance="cov.nug", theta=1))

tau <- 0.499
nx <- 4
ny <- nx
field <- rgrf(1, nx=nx, ny=ny, tau=tau, Covariance="cov.nug", theta=c(.2, 1, 1.5))
plot(attr(field,"locs"))
quilt.plot(cbind(attr(field,"locs"),z=field), nx=nx, ny=ny)
points(attr(field,"locs"))

# unit grid:
nx <- 12
ny <- 7
field <- rgrf(1, nx=nx, ny=ny, xlim=c(0,nx), ylim=c(0,ny),  Covariance="cov.nug", theta=1)
quilt.plot(cbind(attr(field,"locs"),z=field), nx=nx, ny=ny)
points(attr(field,"locs"))

# effective tau:
tau <- 0.25
nx <- 12
ny <- 7
field <- rgrf(1, nx=nx, ny=ny, xlim=c(0,nx), ylim=c(0,ny), tau=tau,  Covariance="cov.nug", theta=1)
quilt.plot(cbind(attr(field,"locs"),z=field), nx=nx, ny=ny)
points(attr(field,"locs"))
attr(field,"eff")

tau <- 0.25
nx <- 12
ny <- 7
field <- rgrf(1, nx=nx, ny=ny, xlim=c(0,nx), ylim=c(0,ny), tau=tau,  Covariance="cov.nug", theta=1)
quilt.plot(cbind(round(attr(field,"locs")+.499)-.49,z=field), nx=nx, ny=ny)
points(attr(field,"locs"))
attr(field,"eff")

# min distance checking:  (within 0,1)
tau <- 0.499
nx <- 40
ny <- nx
field <- rgrf(1, nx=nx, ny=ny, tau=tau, Covariance="cov.nug", theta=c(.2, 1, 1.5))

tt <- rdist(attr(field,"locs"))
etau <- attr(field,"effective.tau")[1]

expect_lt( (1/nx-2*etau), min( tt[tt>0]))
expect_lt( max( tt[tt>0]), sqrt(2)*((nx-1)/nx+2*etau))



context("test-rgrf.R: mean")
tau <- 0.499
nx <- 20
ny <- nx
field <- rgrf(1, nx=nx, ny=ny, tau=tau, Covariance="cov.nug", theta=c(.2, 1, 1.5))

locs <- attr(field,"locs")
field <- rgrf(1, nx=nx, ny=ny, tau=tau, X=cbind(1,locs), beta=c(4,1,1), Covariance="cov.nug", theta=c(.2, 1, 1.5))

quilt.plot(cbind(attr(field,"locs"),z=field))



context("test-rgrf.R: samples")


nx <- 12
ny <- nx
set.seed(12)
field <- rgrf(1, nx=nx, ny=ny, Covariance="cov.sph", theta=c(1,.2,.2))

distmat <- rdist( attr(field,"locs"))

set.seed(12)
out <- rmvnorm(1, Sigma=cov.sph(distmat,theta=c(1,.2,.2)))

expect_equal( norm(field-out), 0)


expect_error( rgrf(1, nx=nx, ny=ny, Covariance="cov.sphd", theta=c(1,.2,.2)))
expect_error( rgrf(1, nx=nx, ny=ny, Covariance="cov.sphd", theta=c(1,.2,.2)))
expect_error( rgrf(1, nx=nx, ny=ny, Covariance="cov.sphd", theta=c(1,.2,.2)))
expect_error( rgrf(1, nx=nx, ny=ny, Covariance="cov.sphd", theta=c(1,.2,.2)))




if (F) {
    n <- 2000
    locs <- cbind(runif(n), runif(n))

    # rdist just beats the shit out of it...
#    microbenchmark::microbenchmark(
#                        a=rdist(locs),
#                        b=as.matrix(dist(locs)),
#                        c=nearest.dist(locs, upper=NULL, delta=.5),
#                        d=nearest.dist(locs, upper=NULL, delta=sqrt(2)),
#                        times=10)

    }
