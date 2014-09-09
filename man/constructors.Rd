% This is file ../spam/man/constructors.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   

\name{constructors}
\alias{rowpointers}
\alias{rowpointers<-}
\alias{colindices}
\alias{colindices<-}
\alias{entries}
\alias{entries<-}
\alias{dimension<-}
\alias{constructors}
\title{Slot modification}
\description{Modify slots of \code{spam} objects
}
\usage{
rowpointers( x) <- value
colindices( x) <- value
entries( x) <- value}
\arguments{
  \item{x}{a \code{spam} matrix}
  \item{value}{vector of appropriate length.}
  }
\value{Modified \code{spam} object.}
\details{Various tests are performed. Thus much slower than direct
  assignment.\cr
  Slot \code{dimension} should be changed through \code{pad} or \code{dim}
}
  \examples{
x <- diag.spam( 2)  
rowpointers( x) <- c(1,1,3)

# The last line is equivalent to 
x@rowpointers <- as.integer( c(1,1,3))
}
\author{Reinhard Furrer}
\keyword{array}
