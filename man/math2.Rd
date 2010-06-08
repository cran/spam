% This is file ../spam0.22-0/man/math2.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{Math2}
\alias{Math2}
\alias{Math2,spam-method}
\alias{Math2,spam,numeric-method}
\alias{round}
\alias{signif}
\title{Rounding of Numbers}
\description{Applies the \code{Math2} group functions to '\code{spam}' objects
}
\usage{# round(x, digits = 0)
# signif(x, digits = 6)
}
\arguments{\item{x}{spam object.}
\item{digits}{integer indicating the precision to be used.}
}
\value{All functions operate on the vector \code{x@entries} and return the
  result thereof.
}
\details{For this generic class typical \code{na.rm} argument has
    no weight here as NA/NaN/Inf are not meaningful (yet).\cr


%  Is implemented for R>=2.3.x only. However,
%  it would be possible to use Martin's proposed workaround:
%    \url{http://tolstoy.newcastle.edu.au/R/help/05/12/18192.html}
}
%\references{
%}
\seealso{\code{Math}}
\examples{
getGroupMembers("Math2")

smat <- diag.spam( rnorm(15))
round(smat, 3)
}
\author{Reinhard Furrer}
\keyword{manip}
