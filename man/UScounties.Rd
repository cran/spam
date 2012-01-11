% This is file ../spam0.28-0/man/UScounties.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{UScounties}
\alias{UScounties}
\alias{UScounties.storder}
\alias{UScounties.ndorder}
\docType{data}
\title{
Adjacency structure of the counties in the contiguous United States 
}
\description{
  First and second order adjacency structure of the counties
  in the contiguous United States. We consider that two counties
  are neighbors if they share at
  least one edge of their polygon description in \code{maps}. 
}

\format{
  Two matrices of class \code{spam}

\describe{
\item{UScounties.storder}{
  Contains a one in the \code{i} and \code{j} element if  county
  \code{i} is a neighbor of county \code{j}.
}
\item{UScounties.ndorder}{
  Contains a one in the \code{i} and \code{j} element if  counties
  \code{i} and \code{j} are a neighbors of county \code{k} and counties 
  \code{i} and \code{j}  are not  neighbors.
}
}
}
%\source{\url{ www.to somethin??}}
\seealso{\code{map} from \pkg{maps}.}
%\references{}
\examples{
# number of counties:
n  <- nrow( UScounties.storder)

\dontrun{
# make a precision matrix 
Q <- diag.spam( n) + .2 * UScounties.storder + .1 * UScounties.ndorder
display( as.spam( chol( Q)))
}
} 
\keyword{datasets}
