% This is file ../spam0.21-0/man/rmvnorm.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{rmvnorm}
\alias{rmvnorm.spam}
\alias{rmvnorm.prec}
\alias{rmvnorm.canonical}
\title{Draw Multivariate Normals}
\description{
  Fast ways to draw multivariate normals when the variance or precision matrix
  is sparse.}
\usage{
rmvnorm.spam(n,mu=rep(0, nrow(Sigma)), Sigma, Rstruct=NULL, ...) 
rmvnorm.prec(n,mu=rep(0, nrow(Q)), Q, Rstruct=NULL, ...) 
rmvnorm.canonical(n, b, Q, Rstruct=NULL, ...) 
}
\arguments{
\item{n}{number of observations.}
\item{mu}{mean vector.}
\item{Sigma}{covariance matrix of class \code{spam}.}
\item{Q}{precision matrix.}
\item{b}{vector determining the mean.}
\item{Rstruct}{the Cholesky structure of \code{Sigma} or \code{Q}.}
\item{\dots}{arguments passed to \code{chol}.}
}
\details{The functions \code{rmvnorm.prec} and \code{rmvnorm.canonical}
  do not requrie sparse precision matrices. 
  For \code{rmvnorm.spam}, the differences between regular and sparse
  covariance matrices are too significant to be implemented here.  \cr
  Often (e.g., in a Gibbs sampler setting), the sparsity structure of
  the covariance/precision does not change. In such setting, the
  Cholesky factor can be passed via \code{Rstruct} in which only updates
  are performed (i.e., \code{update.spam.chol.NgPeyton} instead of a
  full \code{chol}).
}
%\note{There is intentionally no \acronym{S3} distinction between the classes
%  \code{spam} and \code{spam.chol.}\emph{method}.}
\references{See references in \code{\link{chol}}. 
}
\seealso{\code{\link{chol}} and \code{\link{ordering}}.
}
\examples{
# Generate multivariate form a covariance inverse:
# (usefull for GRMF)
set.seed(13)
n <- 25    # dimension
N <- 1000  # sample size
Sigmainv <- .25^abs(outer(1:n,1:n,"-"))
Sigmainv <- as.spam( Sigmainv, eps=1e-4)


Sigma <- solve( Sigmainv)  # for verification 
iidsample <- array(rnorm(N*n),c(n,N))

mvsample <- backsolve( chol(Sigmainv), iidsample)
norm( var(t(mvsample)) - Sigma, type="HS")

# compare with:
mvsample <- backsolve( chol(as.matrix( Sigmainv)), iidsample)
norm( var(t(mvsample)) - Sigma, type="HS")



# 'solve' step by step:
b <- rnorm( n)
R <- chol(Sigmainv)
norm( backsolve( R, forwardsolve( R, b))-
      solve( Sigmainv, b),type="HS") 
norm( backsolve( R, forwardsolve( R, diag(n)))- Sigma,type="HS") 


}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer}
\keyword{algebra}
