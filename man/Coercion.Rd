% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/coercion.Rd.                        %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\name{Coercion}
\docType{class}
\alias{as.vector}
\alias{as.vector.spam}
\alias{as.vector,spam-method}
\alias{as.vector,spam.chol.NgPeyton-method}
\title{Coercion to a Vector}
\description{Coercion of \code{spam} matrices to proper vector objects
}
\usage{\S4method{as.vector}{spam}(x, mode = "any")
}
\arguments{\item{x}{spam object.}
\item{mode}{character string naming an atomic mode or \code{"any"}/\code{"list"}/\code{"expression"}.}
}
\value{If \code{structurebased=TRUE}, the vector \code{x@entries}.\cr
  Conversely, if \code{structurebased=FALSE}, the result is identical to
  one with \code{as.vector(as.matrix(x))}.
}
\details{This coercion allows smooth transitions between different
  matrix formats, see example below.\cr
The Cholesky factors are first transformed to a \code{spam} object.
}
%\references{
%}
\seealso{\code{\link{spam.options}}}
\examples{
x <- diag(2)
ifelse( x, x, 1-x)
ifelse( x, as.vector(x), 1-as.vector(x))

x <- diag.spam(2)
options(spam.structurebased=FALSE)
ifelse( x, as.vector(x), 1-as.vector(x))
options(spam.structurebased=TRUE)
ifelse( x, as.vector(x), 1-as.vector(x))
}
\author{Reinhard Furrer}
\keyword{manip}
