% This is file ../spam0.29-1/man/math.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{Math}
\alias{Math}
\alias{Math,spam-method}
\alias{ceiling}
\alias{floor}
\alias{trunc}
\alias{exp}
\alias{log}
\alias{log10}
\alias{sqrt}
\alias{gamma}
\alias{lgamma}
\alias{abs}
\alias{cumprod}
\alias{cumsum}
\alias{cos}
\alias{acos}
\alias{acosh}
\alias{sin}
\alias{asin}
\alias{asinh}
\alias{tan}
\alias{atan}
\alias{atanh}
\title{Mathematical functions}
\description{Applies the \code{Math} group functions to '\code{spam}' objects
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
...
}
\arguments{\item{x}{spam object.}
\item{base}{positive number.  The base with respect to which logarithms
          are computed.  Defaults to \code{e=exp(1)}.}
}
\value{All functions operate on the vector \code{x@entries} and return the
  result thereof.
}
\details{
  It is important to note that the zero entries do not enter the
  evaluation. The operations are performed on the stored non-zero
  elements. This may lead to differences if compared with the same
  operation on a full matrix. For example, the cosine of sparse matrix is a
  (full) matrix with many ones.

  Evaluating function resulting in NA/NaN/Inf is possible but the result
  cannot be used further as NA/NaN/Inf are not meaningful (yet).\cr
}
%\references{
%}
\seealso{\code{\link{Math2}}
}
\examples{
getGroupMembers("Math")

mat <- matrix(c( 1,2,0,3,0,0,0,4,5),3)
smat <- as.spam( mat)
cos( mat)
cos( smat)


sqrt( smat)
}
\author{Reinhard Furrer}
\keyword{manip}
