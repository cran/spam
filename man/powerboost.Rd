% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/powerboost.Rd.                      %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      %
%    Roman Flury [ctb], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



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

