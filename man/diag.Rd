% This is file ../spam/man/diag.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   










\name{diag}
\alias{diag}
\alias{diag.of.spam}
\alias{diag.spam}
\alias{diag,spam-method}
\alias{diag,ANY-method} 
\alias{diag<-} 
\alias{diag<-,ANY-method} 
\alias{diag<-,spam-method}
\alias{diag<-.spam}
\alias{diag.spam<-}
\alias{diag.assign,spam-method} 
\title{Sparse Matrix diagonals}
\description{Extract or replace the diagonal of a matrix, or construct a
     diagonal matrix.
}
\usage{
# diag(x)
diag(x=1, nrow, ncol)
diag(x) <- value

diag.spam(x=1, nrow, ncol)
diag.spam(x) <- value
}
\arguments{
  \item{x}{a \code{spam} matrix, a vector or a scalar.}
\item{nrow, ncol}{Optional dimensions for the result.}
\item{value}{either a single value or a vector of length equal to that of
   the current diagonal.}
 }
\value{
   If \code{x} is a spam matrix then \code{diag(x)} returns the diagonal of \code{x}.

   The assignment form sets the diagonal of the sparse matrix \code{x} to the
   given value(s).\cr

   \code{diag.spam} works as \code{diag} for spam matrices:
   If \code{x} is a vector (or 1D array) of length two or more, then
   \code{diag.spam(x)} returns a diagonal matrix whose diagonal is \code{x}.

   If \code{x} is a vector of length one then \code{diag.spam(x)} returns an
   identity matrix of order the nearest integer to \code{x}.  The
   dimension of the returned matrix can be specified by \code{nrow} and
   \code{ncol} (the default is square).

   The assignment form sets the diagonal of the matrix \code{x} to the
   given value(s).

  
}
\details{Using \code{diag(x)} can have unexpected effects if \code{x} is a vector
     that could be of length one. Use \code{diag(x, nrow = length(x))} for
     consistent behaviour.

}
%\references{}
\seealso{\code{\link{upper.tri}}, \code{\link{lower.tri}}.
}
\examples{
diag.spam(2, 4)           # 2*I4
smat <- diag.spam(1:5)
diag( smat)
diag( smat) <- 5:1

# The last line is equivalent to 
diag.spam( smat) <- 5:1

# Note that diag.spam( 1:5) <- 5:1 not work of course.


}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
