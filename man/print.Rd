% This is file ../spam0.15-4/man/print.Rd
% This file is part of the spam package, 
%      http://www.mines.edu/~rfurrer/software/spam/
% written and maintained by Reinhard Furrer.






\name{print}
\alias{print}
\alias{print.spam}
\alias{print,spam-method}
\alias{print,spam.chol.NgPeyton-method}
\alias{print.spam.chol.NgPeyton}
\alias{summary}
\alias{summary.spam}
\alias{summary,spam-method}
\alias{summary,spam.chol.NgPeyton-method}
\alias{summary.spam.chol.NgPeyton}
\title{Printing and summarizing sparse matrices}
\description{Printing (non-zero elements) of sparse matrices
  and summarizing the sparseness structure thereof. 
}
\usage{
  print(x, ...)
  summary(object, ...)
}
\arguments{\item{x}{matrix of class \code{spam} or \code{spam.chol.}\emph{method}.}
  \item{object}{matrix of class \code{spam} or \code{spam.chol.}\emph{method}.}
  \item{...}{any other arguments passed
    to \code{print.default}.}
}
\value{\code{NULL} for \code{print}, because the information is printed with \code{cat} there
  is no real need to pass any object back. \cr
  %
  A list containing the non-zero elements and the density for
  \code{summary}
  for class \code{spam}.\cr
  %
  A list containing the non-zero elements of the factor, the density and
  the fill-in for  \code{summary} for class \code{spam.chol.NgPeyton}.}
\details{\code{spam.getOption('printsize')} determines  if the sparse matrix is
  coerced into a matrix and the printed as an array or if
  only the non-zero elements of the matrix are given.
}
%\references{}
\seealso{\code{\link{spam.options}}}
\examples{
set.seed(13)
nz <- 8
ln <- nz
smat <- spam(0,ln,ln)
smat[cbind(sample(ln,nz),sample(ln,nz))] <- 1:nz

par(mfcol=c(1,2),pty='s')
spam.options( printsize=1000)
print(smat)
spam.options( printsize=10)
print(smat)
summary(smat)
(summary(smat))
}
\author{Reinhard Furrer}
\keyword{hplot}
