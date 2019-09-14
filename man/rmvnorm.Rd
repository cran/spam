% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/rmvnorm.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{rmvnorm}
\alias{rmvnorm.spam}
\alias{rmvnorm.prec}
\alias{rmvnorm.canonical}
\title{Draw Multivariate Normals}
\description{
  Fast ways to draw multivariate normals when the variance or precision matrix
  is sparse.}
\usage{
rmvnorm.spam(n,mu=rep.int(0, dim(Sigma)[1]), Sigma, Rstruct=NULL, ...)
rmvnorm.prec(n,mu=rep.int(0, dim(Q)[1]), Q, Rstruct=NULL, ...)
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
  do not require sparse precision matrices.
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
# Generate multivariate from a covariance inverse:
# (usefull for GRMF)
set.seed(13)
n <- 25    # dimension
N <- 1000  # sample size
Sigmainv <- .25^abs(outer(1:n,1:n,"-"))
Sigmainv <- as.spam( Sigmainv, eps=1e-4)


Sigma <- solve( Sigmainv)  # for verification
iidsample <- array(rnorm(N*n),c(n,N))

mvsample <- backsolve( chol(Sigmainv), iidsample)
norm( var(t(mvsample)) - Sigma, type="m")

# compare with:
mvsample <- backsolve( chol(as.matrix( Sigmainv)), iidsample, n)
   #### ,n as patch
norm( var(t(mvsample)) - Sigma, type="m")



# 'solve' step by step:
b <- rnorm( n)
R <- chol(Sigmainv)
norm( backsolve( R, forwardsolve( R, b))-
      solve( Sigmainv, b) )
norm( backsolve( R, forwardsolve( R, diag(n)))- Sigma )


}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer}
\keyword{algebra}
