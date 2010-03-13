% This is file ../spam0.21-0/man/isSymmetric.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   


\name{isSymmetric}
\alias{isSymmetric.spam}
\alias{isSymmetric,spam-method}
\title{Test if a spam matrix is Symmetric}
\description{Efficient function to test if 'object' is symmetric or not.}
\usage{
isSymmetric.spam(object, tol = 100 * .Machine$double.eps, ...)}
\arguments{
  \item{object}{a \code{spam} matrix.}
\item{tol}{numeric scalar >= 0.  Smaller differences are not considered,
          see \code{all.equal.spam}.}
\item{...}{further arguments passed to \code{all.equal.spam}.}
 }
 \details{symmetry is assessed by comparing the sparsity structure of
   \code{object} and \code{t(object)} via the function
   \code{all.equal.spam}.}
 \value{
   logical indicating if \code{object} is symmetric or not.}
 \seealso{\code{\link{all.equal.spam}}.}
 \examples{
obj <- diag.spam(2)
isSymmetric(obj)

obj[1,2] <- .Machine$double.eps
isSymmetric(obj)
}
 \author{Reinhard Furrer}
\keyword{array}
