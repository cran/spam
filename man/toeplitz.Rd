% This is file ../spam0.28-0/man/toeplitz.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   



\name{toeplitz}
\alias{toeplitz.spam}
\title{Create Toeplitz Matrices}
\description{Creates symmetric and asymmetric Toeplitz matrices.}
\usage{
toeplitz.spam(x, y = NULL, eps = .Spam$eps)
}
\arguments{
  \item{x}{the first row to form the Toeplitz matrix.}
  \item{y}{for asymmetric Toeplitz matrices, this contains the first
    column.}
  \item{eps}{A tolerance parameter:  elements of \code{x} such that \code{abs(x)
    <= eps} set to zero. Defaults to \code{eps = .Spam$eps}.} 
}
\value{The Toeplitz matrix in \code{spam} format.}
\details{The vector \code{y} has to be of the same length as \code{x}
  and its first element is discarded. 
    }
%\references{}
\seealso{\code{\link{toeplitz}}, \code{\link{circulant.spam}}}
\examples{
toeplitz.spam(c(1,.25,0,0,0))

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
