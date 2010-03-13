% This is file ../spam0.21-0/man/version.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{version}
\alias{version}
\alias{spam.version}
\alias{spam.Version}
\title{Spam Version Information}
\description{
  \code{spam.version} is a variable (\code{list}) holding
  detailed information about the version of
  \code{spam} loaded.

  \code{spam.Version()} provides detailed information about the version of \code{spam}
     running.

  
}
\usage{
spam.version
}
\value{\code{spam.version} is a list with character-string components
  \item{status}{the status of the version (e.g., \code{"beta"})}

  \item{major}{the major version number}
  
  \item{minor}{the minor version number}
  
  \item{year}{the year the version was released}
  
  \item{month}{the month the version was released}
  
  \item{day}{the day the version was released}
  
  \item{version.string}{a \code{character} string concatenating
     the info  above, useful for plotting, etc.}

  \code{spam.version} is a list of class \code{"simple.list"} which
     has a \code{print} method.
}
% \references{}
\seealso{See the R counterparts  \code{\link[base]{R.version}}.
}
\author{Reinhard Furrer}
\examples{
spam.version$version.string
}
\keyword{environment}
\keyword{sysdata}
\keyword{programming}
