% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/s3only.Rd.                          %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

\name{s3only}
\alias{var.spam}
\title{Wappers for Sparse Matrices}
\description{
  These functions are convenient wrappers for  \code{spam} objects
  to classical matrix operations.}
\usage{
var.spam(x, \dots)

\S3method{var}{spam}(x, \dots)
}
\arguments{
\item{x}{matrix of class \code{spam}.}
\item{\dots}{further arguments passed to or from other methods.}
}
\value{Depends on function\dots
}
\details{There is probably no point in fully defining methods
  here. Typically, these functions do not exploit sparsity
  structures. Hence, for very large matrices, warnings may be posted.}
%\references{}
%\note{}
\seealso{Option \code{"inefficiencywarning"} in \code{\link{spam.options}}.
}
\author{Reinhard Furrer}
\keyword{algebra}
