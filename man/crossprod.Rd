% This is file ../spam/man/crossprod.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   



\name{crossprod}
\alias{crossprod.default}
\alias{crossprod.spam}
\alias{tcrossprod.spam}
\title{Spam Matrix Crossproduct}
\description{
     Given matrices \code{x} and \code{y} as arguments, return a matrix
     cross-product.  This is formally equivalent to (but usually
     slightly faster than) the call \code{t(x) \%*\% y} (\code{crossprod.spam}) or \code{x
     \%*\% t(y)} (\code{tcrossprod.spam}).}

\usage{

     crossprod.spam(x, y = NULL)
     
     tcrossprod.spam(x, y = NULL)
}
\arguments{
  \item{x, y}{matrices: \code{y = NULL} is taken to be the
          same matrix as \code{x}.  Vectors are promoted to single-column or
          single-row matrices, depending on the context.}
}
\value{A double matrix}
\note{
     When \code{x} or \code{y} are not matrices, they are treated as column or
     row matrices.
   }
%\references{}
%\seealso{\code{\link{chol}}}
\examples{
crossprod.spam(diag.spam(2),1:2)

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
