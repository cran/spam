% This is file ../spam0.15-4/man/summary.Rd
% This file is part of the spam package, 
%      http://www.mines.edu/~rfurrer/software/spam/
% written and maintained by Reinhard Furrer.






\name{Summary}
\alias{Summary}
\alias{Summary.spam}
\alias{Summary,spam-method}
\alias{all}
\alias{any}
\alias{max}
\alias{min}
\alias{prod}
\alias{range}
\alias{sum}
\title{Rounding of Numbers}
\description{Applies the \code{Math2} group functions to '\code{spam}' objects
}
\usage{# max(x,..., na.rm = FALSE)
}
\arguments{\item{x}{spam object.}
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
\seealso{\code{\link{Math}} and \code{\link{Math2}}.}
\examples{
getGroupMembers("Summary")

smat <- diag.spam( runif(15))
range(smat)
}
\author{Reinhard Furrer}
\keyword{manip}
% all any max min prod range sum
