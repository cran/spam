% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/summary.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      %
%    Roman Flury [ctb], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


\name{Summary}
\alias{Summary.spam}
\alias{Summary,spam-method}
\alias{all.spam}
\alias{any.spam}
\alias{max.spam}
\alias{min.spam}
\alias{prod.spam}
\alias{range.spam}
\alias{sum.spam}
\alias{all,spam-method}
\alias{any,spam-method}
\alias{max,spam-method}
\alias{min,spam-method}
\alias{prod,spam-method}
\alias{range,spam-method}
\alias{sum,spam-method}
\title{Rounding of Numbers}
\description{Applies the \code{Math2} group functions to \code{spam} objects
}
\usage{# max(x,..., na.rm = FALSE)
}
\arguments{\item{x}{spam object.}
   \item{na.rm}{a logical indicating whether missing values should be
          removed.}
}
\value{If \code{structurebased=TRUE}, all functions operate on the vector \code{x@entries} and return the
  result thereof.\cr
  Conversely, if \code{structurebased=FALSE}, the result is identical to
  one with \code{as.matrix(x)} input.
}
\details{The \code{na.rm} argument is only meaninful if \code{NAOK=TRUE}.
}
%\references{
%}
\seealso{\code{\link{Math.spam}} and \code{\link{Math2}}.}
\examples{
getGroupMembers("Summary")

smat <- diag.spam( runif(15))
range(smat)
options(spam.structurebased=FALSE)
range(smat)

\dontrun{
max( log(spam(c(1,-1))), na.rm=TRUE)
}
# allow 'NA's first:
# TODO
# options(spam.NAOK=TRUE)
# max( log(spam(c(1,-1))), na.rm=TRUE)
}
\author{Reinhard Furrer}
\keyword{manip}
% all any max min prod range sum
