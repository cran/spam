% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/math2.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


\name{Math2}
\alias{Math2.spam}
\alias{Math2,spam-method}
\alias{Math2,spam,numeric-method}
\alias{round,spam-method}
\alias{signif,spam-method}
\alias{round.spam}
\alias{signif.spam}
\title{Rounding of Numbers}
\description{Applies the \code{Math2} group functions to '\code{spam}' objects
}
\usage{\S4method{round}{spam}(x, digits = 0)
\S4method{signif}{spam}(x, digits = 6)
}
\arguments{\item{x}{spam object.}
\item{digits}{integer indicating the precision to be used.}
}
\value{All functions operate on the vector \code{x@entries} and return the
  result thereof.
}
%\details{%  Is implemented for R>=2.3.x only. However,
%  it would be possible to use Martin's proposed workaround:
%    \url{http://tolstoy.newcastle.edu.au/R/help/05/12/18192.html}
%\references{
%}
\seealso{\code{\link{Ops.spam}} and \code{\link{Math.spam}}}
\examples{
getGroupMembers("Math2")

smat <- diag.spam( rnorm(15))
round(smat, 3)
}
\author{Reinhard Furrer}
\keyword{manip}
% "round"  "signif"
