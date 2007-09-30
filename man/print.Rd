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
\value{\code{NULL}, because the information is printed with \code{cat} there
  is no real need to pass any object back. We also do not see the need
  of defining a \code{summary.spam} class with a proper \code{print}
  method.}
\details{\code{.Spam$printsize} determines  if the sparse matrix is
  coerced into a matrix and the printed as an array or if
  only the non-zero elements of the matrix are given.
}
%\references{}
\seealso{\code{\link{.Spam}}}
}
\examples{
nz <- 8
ln <- nz
smat <- spam(0,ln,ln)
smat[cbind(sample(ln,nz),sample(ln,nz))] <- 1:nz

par(mfcol=c(1,2),pty='s')
.Spam$printsize <- 1000
print(smat)
.Spam$printsize <- 10
print(smat)
summary(smat)

}
\author{Reinhard Furrer}
\keyword{hplot}
