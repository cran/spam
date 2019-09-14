% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/diff.Rd.                            %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

\name{diff}
\alias{diff.spam}
\alias{diff,spam-method}
\title{Lagged Differences}
\description{Returns suitably lagged and iterated differences.}
\usage{
# diff.spam(x, lag = 1, differences = 1, ...)
\S4method{diff}{spam}(x, lag = 1, differences = 1, ...)
}
\arguments{
  \item{x}{a \code{spam} matrix containing the values to be
    differenced.}
  \item{lag}{an integer indicating which lag to use.}
  \item{differences}{an integer indicating the order of the difference.}
  \item{...}{further arguments to be passed to or from methods.}
}
\value{A \code{spam} matrix with elements similar to
  \code{as.spam(diff(as.matrix(x), ...))}.
}
%\details{    }
%\references{}
\seealso{\code{diff} in \code{base}, \code{\link{precmat}}.}
\examples{
# incidence matrix for a RW(3) model
D <- diff.spam(diag.spam(10), lag=1, differences=3)
 t(D)%*%D
}
\author{Reinhard Furrer}
\keyword{array}
\keyword{manip}
