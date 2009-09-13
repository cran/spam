% This is file ../spam0.15-5/man/adiag.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.








\name{adiag}
\alias{adiag}
\alias{adiag.spam}
\title{Binds Arrays Corner-to-Corner}
\description{Creates a sparse block-diagonal matrix.
}
\usage{
adiag.spam(...)
}
\arguments{
  \item{...}{Arrays to be binded together}
 }

\details{
This is a small helper function to create block diagonal sparse matrices. In the two matrix case,
\code{adiag.spam(A,B)}, this is equivalent to a complicated \code{rbind(cbind(A, null), cbind(B, t(null)))},
where \code{null} is a null matrix of appropriate dimension.\cr

It is recursively defined.

The arrays are coerced to sparse matrices first.\cr


This function is similar to the function \code{adiag} from the package \code{magic}.
	      However, here no padding is done and all the dimnames are stripped.
}
\value{
    Returns a \code{spam} matrix as described above.
}
%\references{}
\seealso{
	\code{\link{diag.spam}}.
}
\examples{
A <- diag.spam(2, 4)           # 2*I4
B <- matrix(1,3,3)
AB <- adiag.spam(A,B)

# equivalent to:
ABalt <- rbind(cbind( A, matrix(0,nrow(A),ncol(B))),
               cbind( matrix(0,nrow(B),ncol(A)), B))
         
norm(AB-ABalt)


# Matrices do not need to be square:
adiag.spam(1,2:5,6)
}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
