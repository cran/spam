% This is file ../spam0.29-1/man/spam.chol.NgPeyton-class.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{spam.chol.NgPeyton-class}
\docType{class}
\alias{spam.chol.NgPeyton-class}
\alias{as.matrix,spam.chol.NgPeyton-method}
\alias{as.spam,spam.chol.NgPeyton-method}
\alias{backsolve,spam.chol.NgPeyton-method}
\alias{c,spam.chol.NgPeyton-method}
\alias{determinant,spam.chol.NgPeyton-method}
\alias{diag,spam.chol.NgPeyton-method}
\alias{dim<-,spam.chol.NgPeyton-method}
\alias{dim,spam.chol.NgPeyton-method}
\alias{display,spam.chol.NgPeyton-method}
\alias{forwardsolve,spam.chol.NgPeyton-method}
\alias{image,spam.chol.NgPeyton-method}
\alias{length<-,spam.chol.NgPeyton-method}
\alias{length,spam.chol.NgPeyton-method}
\alias{ordering,spam.chol.NgPeyton-method}
\alias{print,spam.chol.NgPeyton-method}
\alias{show,spam.chol.NgPeyton-method}
\alias{summary,spam.chol.NgPeyton-method}
\alias{t,spam.chol.NgPeyton-method}
\alias{chol,spam.chol.NgPeyton-method}

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
    \item{as.matrix}{\code{signature(x = "spam.chol.NgPeyton")}: Transform the factor
      into a regular matrix. }
    \item{as.spam}{\code{signature(x = "spam.chol.NgPeyton")}:  Transform the factor
      into a \code{spam} object.}
    \item{backsolve}{\code{signature(r = "spam.chol.NgPeyton")}: solving
    a triangular system, see \code{\link{solve}}. }
    \item{forwardsolve}{\code{signature(l = "spam.chol.NgPeyton")}: solving
    a triangular system, see \code{\link{solve}}. }
    \item{c}{\code{signature(x = "spam.chol.NgPeyton")}: Coerce the factor into a vector. }
    \item{determinant}{\code{signature(x = "spam.chol.NgPeyton")}:
      Calculates the determinant from the factor, see also \code{\link{det}}. }
    \item{diag}{\code{signature(x = "spam.chol.NgPeyton")}: Extracts the
    diagonal entries.}
%    \item{dim<-}{\code{signature(x = "spam.chol.NgPeyton")}: ... }
    \item{dim}{\code{signature(x = "spam.chol.NgPeyton")}: Retrieve the
      dimension. Note that \code{"dim<-"} is not implemented.}
    \item{display}{\code{signature(x = "spam.chol.NgPeyton")}: Transformation
      to a \code{spam} object and display, see also \code{\link{display}}. }
    \item{image}{\code{signature(x = "spam.chol.NgPeyton")}: Transformation
      to a \code{spam} object and display, see also \code{\link{image}}. }
%    \item{length<-}{\code{signature(x = "spam.chol.NgPeyton")}: ... }
    \item{length}{\code{signature(x = "spam.chol.NgPeyton")}: Retrieve the
      dimension. Note that \code{"length<-"} is not implemented. }
    \item{ordering}{\code{signature(x = "spam.chol.NgPeyton")}:
      Retrieves the ordering, in \code{\link{ordering}}. }
    \item{print}{\code{signature(x = "spam.chol.NgPeyton")}: Short description. }
    \item{show}{\code{signature(object = "spam.chol.NgPeyton")}: Short description. }
    \item{summary}{\code{signature(object = "spam.chol.NgPeyton")}:
      Description of the factor, returns (as a list) \code{nnzR}, \code{nnzcolindices},
      the density of the factor \code{density}, and fill-in ratio
      \code{fillin}. For the use of the first two, see \sQuote{Examples}
    in \code{\link{chol}}.}
    \item{t}{\code{signature(x = "spam.chol.NgPeyton")}: Transformation
      to a \code{spam} object and transposition. }
    \item{chol}{\code{signature(x = "spam.chol.NgPeyton")}: Returns
      \code{x} unchanged. }
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
