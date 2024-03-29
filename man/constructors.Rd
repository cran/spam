% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/constructors.Rd.                    %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

\name{constructors}
\alias{rowpointers}
\alias{rowpointers<-}
\alias{colindices}
\alias{colindices<-}
\alias{entries}
\alias{entries<-}
\alias{dimension<-}
\alias{constructors}
\title{Slot Modification}
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
