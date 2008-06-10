\name{spam.chol.NgPeyton-class}
\docType{class}
\alias{spam.chol.NgPeyton-class}
\alias{backsolve,spam.chol.NgPeyton-method}
\alias{forwardsolve,spam.chol.NgPeyton-method}
\alias{show,spam.chol.NgPeyton-method}
\alias{dim<-,spam.chol.NgPeyton-method}
\alias{dim,spam.chol.NgPeyton-method}
\alias{diag,spam.chol.NgPeyton-method}
\alias{as.matrix,spam.chol.NgPeyton-method}
\alias{length<-,spam.chol.NgPeyton-method}
\alias{length,spam.chol.NgPeyton-method}
\alias{c,spam.chol.NgPeyton-method}
\title{Class "spam.chol.NgPeyton"}
\description{Result of a Cholesky decomposition with the \code{NgPeyton}
method}
\section{Objects from the Class}{
Objects are created by calls of the form \code{chol(x,method="NgPeyton", ...)}
and should not be created directly with a
\code{new("spam.chol.NgPeyton", ...)} call.\cr
At present, no proper print method is defined. However, the factor can
be transformed into a \code{spam} object.
}
\section{Methods}{
  \describe{
    \item{as.spam}{\code{signature(x = "spam")}: Transform the factor
      into a \code{spam} object  }
    \item{length}{\code{signature(x = "spam")}: ... }
    \item{backsolve}{\code{signature(r = "spam.chol.NgPeyton")}: solving
    a triangular system, see \code{\link{solve}}.}
    \item{forwardsolve}{\code{signature(l = "spam.chol.NgPeyton")}: solving
    a triangular system, see \code{\link{solve}}.}
    \item{dim}{\code{signature(x = "spam")}: Retrieve the
      dimension. Note that \code{"dim<-"} is not implemented.}
    \item{length}{\code{signature(x = "spam")}: Retrieve the
      dimension. Note that \code{"dim<-"} is not implemented. }
    \item{c}{\code{signature(x = "spam")}: Coerce the factor into a vector .}
  }
}
\references{Ng, E. G. and B. W. Peyton (1993), "Block sparse Cholesky algorithms
on advanced uniprocessor computers", \emph{SIAM J. Sci. Comput.}, \bold{14},
pp. 1034-1056.

}
\author{Reinhard Furrer}
\seealso{\code{\link{print.spam}}
  \code{\link{ordering}} and \code{\link{chol}}}
\examples{
x <- spam( c(4,3,0,3,5,1,0,1,4),3)
cf <- chol( x)
cf
as.spam( cf)


# Modify at own risk...
slotNames(cf)

}
\keyword{classes}
