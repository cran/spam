% This is file ../spam0.28-0/man/powerboost.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{powerboost}
\alias{powerboost}
\title{Specific options Setting}
\description{
  Sets several options for speed-up.
}
\usage{
powerboost(flag)
}
\arguments{
  \item{flag}{on or off}
}
\details{The options turn checking off (\code{"safemode"},
  \code{"cholsymmetrycheck"}  and \code{"cholpivotcheck"}) and switch to
  single precision for \code{"eps"}.}
\value{
  \code{NULL} in any case.

}
\seealso{\code{\link{spam.options}}.
}
\author{Reinhard Furrer, after receiving too much C.mc.st adds.}
\keyword{environment}

