% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/pad.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   




\name{pad}
\alias{pad} 
\alias{pad.spam}
\alias{pad<-} 
\alias{pad<-.spam} 
\alias{pad<-,spam-method}
\alias{pad<-,matrix-method}
\title{Padding a (sparse) matrix}
\description{
  Resets the dimension of a (\code{spam}) matrix by truncation or zero padding.
}
\usage{
pad(x) <- value
%"pad<-.spam"(x,value)
}
\arguments{
  \item{x}{a \code{spam} matrix}
  \item{value}{A numeric two-vector.}
}
\value{
  A (\code{spam}) matrix of dimension \code{value} where trunction or
  padding has
  been used.
 
}
\details{
  It is important to notice the different behavior of the replacement
  method for ordinary arrays and \code{spam} objects (see
  \sQuote{Examples}). Here, the elements are not simply rearranged but
  an entirely new matrix is constructed. If the new column dimension is
  smaller than the original, the matrix is also cleaned (with
  \code{option("spam.eps")} as filter).
}
%\references{}
\seealso{\code{\link{dim.spam}}.
}
\examples{
x <- diag(4)
dim(x)<-c(2,8)
x

s <- diag.spam(4)
pad(s) <- c(7,3)  # any positive value can be used

s <- diag.spam(4)
pad(s) <- c(2,8)  # result is different than x
}
\author{Reinhard Furrer}
\keyword{array}
