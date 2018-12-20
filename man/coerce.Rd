% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/coerce.Rd.                          %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      %
%    Roman Flury [ctb], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{coerce-methods}
\docType{methods}
\alias{coerce.spam}
\alias{coerce,spam,logical-method}
\alias{coerce,spam,matrix-method}
\alias{coerce,spam,vector-method}
\alias{coerce,spam,list-method}
\alias{coerce,spam,integer-method}
\title{Force a \code{spam} Object to Belong to a Class}
\description{
  These functions manage the relations that allow coercing a \code{spam} object
  to a given class. 
}
\section{Methods}{
\describe{
\item{\code{signature(from = "spam", to = "matrix")}}{
  this is essentially equivalent to \code{as.matrix(object)}.
}
\item{\code{signature(from = "spam", to = "list")}}{
  this is essentially equivalent to \code{triplet(object)}.
}
\item{\code{signature(from = "spam", to = "vector")}}{
  this is essentially equivalent to \code{object@entries}
  (\code{structurebased=TRUE}) or \code{c(object)}.
}
\item{\code{signature(from = "spam", to = "logical")}}{
  the entries are forced to logicals (nonzeros only in case of
  \code{structurebased=TRUE}).
}
\item{\code{signature(from = "spam", to = "integer")}}{
  the entries are forced to integers (nonzeros only in case of
  \code{structurebased=TRUE}).
}
}}
\examples{
ifelse( diag.spam(2)*c(0,1), TRUE, FALSE)
}
\keyword{methods}