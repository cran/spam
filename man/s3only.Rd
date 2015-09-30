% This is file ../spam/man/s3only.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   

\name{s3only}
\alias{var.spam}
\alias{eigen.spam}
\title{Wappers for Sparse Matrices}
\description{
  These functions are convenient wrappers for  \code{spam} objects
  to classical matrix operations.}
\usage{
var.spam(x, \dots)
eigen.spam(x, \dots)

\S3method{var}{spam}(x, \dots)
\S3method{eigen}{spam}(x, \dots)
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
\examples{
eigen( diag.spam(3))
}
\author{Reinhard Furrer}
\keyword{algebra}