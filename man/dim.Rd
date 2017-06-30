% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/dim.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


\name{dim}
\alias{dim.spam}
\alias{dim<-.spam} 
\alias{dim<-,spam-method}
\title{Dimensions of an Object}
\description{
  Retrieve or set the dimension of an \code{spam} object.
}
\usage{
# dim(x)
# dim(x) <- value
%"dim<-.spam"(x,value)
}
\arguments{
  \item{x}{a \code{spam} matrix}
  \item{value}{A numeric two-vector, which is coerced to integer (by truncation).}
 }
\value{
  \code{dim} retrieves
  the \code{dimension} slot of the object.  It is a vector
  of mode \code{integer}.

  The replacemnt method changes the dimension of the object by rearranging.
  
}
\details{
  In older version of \code{spam}, the behavior of the replacement
  method was different and is now implemented in \code{\link{pad.spam}}.

}
%\references{}
\seealso{\code{\link{pad.spam}}.
}
\examples{
x <- diag(4)
dim(x)<-c(2,8)
x

s <- diag.spam(4)
dim(s) <- c(2,8)  # result is different than x

s <- diag.spam(4)
pad(s) <- c(7,3)  # any positive value can be used

}
\author{Reinhard Furrer}
\keyword{array}
