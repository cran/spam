% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/UScounties.Rd.                      %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
\seealso{\code{map}, from \pkg{maps}.}
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
