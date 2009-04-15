% This is file ../spam0.15-4/man/kronecker.Rd
% This file is part of the spam package, 
%      http://www.mines.edu/~rfurrer/software/spam/
% written and maintained by Reinhard Furrer.






\name{kronecker}
\alias{kronecker}
\alias{kronecker.default}
\alias{kronecker.spam}
\title{Kronecker products on sparse matrices}
\description{Computes the generalised kronecker product of two arrays, \code{X} and
     \code{Y}.}
\usage{
kronecker.spam(X, Y, FUN = "*", make.dimnames = FALSE, ...)
}
\arguments{
  \item{X}{sparse matrix of class \code{spam}, a vector or a matrix.}
  \item{Y}{sparse matrix of class \code{spam}, a vector or a matrix.}
  \item{FUN}{a function; it may be a quoted string. See details}
  \item{make.dimnames}{Provide dimnames that are the product of the dimnames of
    'X' and 'Y'. } 
    \item{ ...}{optional arguments to be passed to \code{FUN}.}
}
\value{An array \code{A} with dimensions
     \code{dim(X) * dim(Y)}.}
\details{The sparseness structure is determined by the ordinary
  \code{\%x\%}. Hence, the result of \code{kronecker(X, Y, FUN = "+")} is
    different depending on the input.
    }
%\references{}
%\seealso{\code{\link{chol}}}
\examples{
# Starting with non-spam objects, we get a spam matrix
kronecker.spam( diag(2), array(1:4,c(2,2)))

kronecker( diag.spam(2), array(1:4,c(2,2)))

# Notice the preservation of sparseness structure: 
kronecker( diag.spam(2), array(1:4,c(2,2)),FUN="+")


}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
