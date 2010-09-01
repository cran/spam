% This is file ../spam0.23-0/man/dim.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{dim}
\alias{dim}
\alias{dim<-.spam} 
\alias{dim<-,spam-method}
\title{Dimensions of an Object}
\description{
  Retrieve or set the dimension of an \code{spam} object.
}
\usage{
# dim(x)
# dim(x) <- value
"dim<-.spam"(x,value)
}
\arguments{
  \item{x}{a \code{spam} matrix}
  \item{value}{A numeric two-vector, which is coerced to integer (by truncation).}
 }
\value{
  \code{dim} retrieves
  the \code{dimension} slot of the object.  It is a vector
  of mode \code{\link{integer}}.

  The replacemnt method changes the dimension of the object by
  truncation or extension (with zeros).
  
}
\details{
  It is important to notice the different behavior of the replacement
  method for ordinary arrays and \code{spam} objects (see
  \sQuote{Examples}). Here, the elements are not simply rearranged but
  an entirely new matrix is constructed. If the new column dimension is
  smaller than the original, the matrix is also cleaned (with
  \code{spam.option("eps")} as filter).

  For the same operation as with regular arrays, use \code{spam}

}
%\references{}
\seealso{\code{\link[base]{dim}}.
}
\examples{
x <- diag(4)
dim(x)<-c(2,8)
x

s <- diag.spam(4)
dim(s) <- c(7,3)  # any positive value can be used

s <- diag.spam(4)
dim(s) <- c(2,8)  # result is different than x

}
\author{Reinhard Furrer}
\keyword{array}
