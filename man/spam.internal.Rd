% This is file ../spam0.15-4/man/spam.internal.Rd
% This file is part of the spam package, 
%      http://www.mines.edu/~rfurrer/software/spam/
% written and maintained by Reinhard Furrer.






\name{spam internal}
\alias{dcheck}
\alias{icheck}
\title{
  Spam internal and auxiliary functions
}
\description{
  The functions listed below are auxiliary functions but are 
  exported by the NAMESPACE. The user should not
  require to call these directly. 
}
\details{The functions are listed here for a better understanding of the
  code (to fulfill the tutorial style paradigm).\cr

  \code{validspamobject(object)} A few sanity checks if \code{object}
  is a proper \code{spam} object.

  \code{dcheck(x)}, \code{icheck(x)} testing and forcing of \code{x}
  to doubles and integers. These functions are used when calling Fortran
  routines.
}
\note{The integers \code{int0}, \code{int1} and \code{int2} (0-2) are not
  exported. 
}
\keyword{internal}
