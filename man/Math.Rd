% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/math.Rd.                            %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




\name{Math}
\alias{Math}
\alias{Math.spam}
\alias{Math,spam-method}
\alias{ceiling.spam}
\alias{ceiling,spam-method}
\alias{floor.spam}
\alias{floor,spam-method}
\alias{trunc.spam}
\alias{trunc,spam-method}
\alias{exp.spam}
\alias{exp,spam-method}
\alias{log.spam}
\alias{log,spam-method}
\alias{log2.spam}
\alias{log2,spam-method}
\alias{log10.spam}
\alias{log10,spam-method}
\alias{sqrt.spam}
\alias{sqrt,spam-method}
\alias{gamma.spam}
\alias{gamma,spam-method}
\alias{digamma.spam}
\alias{digamma,spam-method}
\alias{trigamma.spam}
\alias{trigamma,spam-method}
\alias{lgamma.spam}
\alias{lgamma,spam-method}
\alias{abs.spam}
\alias{abs,spam-method}
\alias{cumprod.spam}
\alias{cumprod,spam-method}
\alias{cumsum.spam}
\alias{cumsum,spam-method}
\alias{cummax.spam}
\alias{cummax,spam-method}
\alias{cummin.spam}
\alias{cummin,spam-method}
\alias{cos.spam}
\alias{cos,spam-method}
%\alias{cospi.spam}
%\alias{cospi,spam-method}
\alias{cosh.spam}
\alias{cosh,spam-method}
\alias{acos.spam}
\alias{acos,spam-method}
\alias{acosh.spam}
\alias{acosh,spam-method}
\alias{sin.spam}
\alias{sin,spam-method}
\alias{asin.spam}
\alias{asin,spam-method}
\alias{asinh.spam}
\alias{asinh,spam-method}
\alias{tan.spam}
\alias{tan,spam-method}
\alias{atan.spam}
\alias{atan,spam-method}
\alias{atanh.spam}
\alias{atanh,spam-method}
\title{Mathematical Functions}
\description{Applies the \code{Math} group functions to \code{spam} objects
}
\usage{# ceiling(x)
# floor(x)
# exp(x, base = exp(1))
# log(x, base = exp(1))
# sqrt(x)

# abs(x)
# cumprod(x)
# cumsum(x)

# cos(x)
# sin(x)
# tan(x)

# acosh(x)

}
\arguments{\item{x}{spam object.}
\item{base}{positive number. The base with respect to which logarithms
are computed. Defaults to \code{e=exp(1)}.}
}
\value{If the option \code{spam.structurebased=TRUE}, all functions operate on the vector \code{x@entries} and return the
result thereof.\cr
Conversely, if \code{structurebased=FALSE}, the result is identical to
one with \code{as.matrix(x)} input and an \code{as.spam} purger.
}
\details{
It is important to note that the zero entries do not enter the
evaluation when \code{structurebased=FALSE}. The operations are performed on the stored non-zero
elements. This may lead to differences if compared with the same
operation on a full matrix.
}
%\references{
%}
\seealso{\code{\link{Summary.spam}}, \code{\link{Ops.spam}} and \code{\link{Math2.spam}}
}
\examples{
getGroupMembers("Math")

mat <- matrix(c( 1,2,0,3,0,0,0,4,5),3)
smat <- as.spam( mat)
cos( mat)
cos( smat)

options(spam.structurebased=FALSE)
cos( smat)

sqrt( smat)
}
\author{Reinhard Furrer}
\keyword{manip}
