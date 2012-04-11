% This is file ../spam0.29-0/man/cleanup.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{cleanup}
\alias{cleanup}
\title{Cleaning up sparse matrices}
\description{
  Eliminates an zeros in a sparse matrix.}
\usage{
cleanup(x, eps = .Spam$eps)
}
\arguments{
\item{x}{a sparse matrix of class \code{spam}.}
\item{eps}{numeric scalar > 0.  Smaller entries are coerced to zero.}
}
\details{A sparse matrix may still contain zeros. This function (aliased
  to \code{as.spam}) filters these values.\cr
  This often causes confusion when testing such matrices for symmetry
  or comparing apparently equal matrices with \code{all.equal}
  (see \sQuote{Examples} below.
}
\seealso{\code{\link{isSymmetric.spam}} and \code{\link{all.equal.spam}}.
}
\examples{
A <- diag.spam(2)
A[1,2] <- 0

all.equal(A, t(A))
isSymmetric.spam(A)
all.equal(cleanup(A), diag.spam(2))
}
\author{Reinhard Furrer}
\keyword{algebra}
