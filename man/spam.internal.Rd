% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/spam.internal.Rd.                   %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   










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
%\note{The integers \code{int0}, \code{int1} and \code{int2} (0-2) are not
%  exported. 
%}
\keyword{internal}
