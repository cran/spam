% This is file ../spam0.29-2/man/precmat.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{precmat}
\alias{precmat}
\alias{precmat.RW1}
\alias{precmat.RW2}
\alias{precmat.season}
\alias{precmat.IGMRFreglat}
\alias{precmat.IGMRFirreglat}
\title{IGMRF Precision Matrices}
\description{
  Fast ways to create sparse precision matrices for various IGMRF.}
\usage{
precmat(n, season=12, m=n, A=NULL, ... , type="RW1")
precmat.RW1(n)
precmat.RW2(n)
precmat.season(n, season=12)
precmat.IGMRFreglat(n, m, anisotropy=1)
precmat.IGMRFirreglat(A, eps=.Spam$eps)
}
\arguments{
\item{n}{dimension of the field.}
\item{type}{the type of the IGMRF.}
\item{season}{length of season.}
\item{m}{second dimension (in case of a regular lattice)}
\item{A}{adjacency matrix (see below)}
\item{anisotropy}{anisotropy factor, between 0 and 2.}
\item{eps}{tolerance level.}
\item{\dots}{arguments passed to individual functions.}
}
\details{\code{precmat} is a wrapper that calls the other functions
  according to the argument \code{type}.
  \cr
  Implements many of the precision matrices discussed in Chapter
  3 of Rue and Held (2005). For example, \code{precmat.RW1},
  \code{precmat.RW2} and \code{precmat.season} are given in
  equations (3.22), (3.40) and (3.59); \code{precmat.IGMRFreglat} on
  page 107. 
  \cr
  If adjacency matrix is a regular matrix, it is coerced to a
  \code{spam} object. Only the structure is used. Make sure, that the
  diagonal is empty.
}
\value{A sparse precision matrix.}
%\note{There is intentionally no \acronym{S3} distinction between the classes
%  \code{spam} and \code{spam.chol.}\emph{method}.}
\references{Rue and Held (2005). 
}
\seealso{\code{\link{rmvnorm.prec}}, \code{\link{adjacency.landkreis}}.
}
\examples{
n <- 10
Q <- precmat.RW2( n)
# rmvnorm.prec(1, Q=Q) # does not work, because the matrix is singular.

Q%*%cbind(1,1:n)
}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer}
\keyword{algebra}