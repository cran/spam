% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/cleanup.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\name{cleanup}
\alias{cleanup}
\title{Cleaning up sparse matrices}
\description{
  Eliminates an zeros in a sparse matrix.}
\usage{
cleanup(x, eps = getOption("spam.eps"))
}
\arguments{
\item{x}{a sparse matrix of class \code{spam}.}
\item{eps}{numeric scalar > 0.  Smaller entries are coerced to zero.}
}
\details{A sparse matrix may still contain zeros. This function (aliased
  to \code{as.spam}) filters these values.\cr
  This often causes confusion when testing such matrices for symmetry
  or comparing apparently equal matrices with \code{all.equal}
  (see \sQuote{Examples} below.
}
\seealso{\code{\link{isSymmetric.spam}} and \code{\link{all.equal.spam}}.
}
\examples{
A <- diag.spam(2)
A[1,2] <- 0

all.equal(A, t(A))
isSymmetric.spam(A)
all.equal(cleanup(A), diag.spam(2))
}
\author{Reinhard Furrer}
\keyword{algebra}
