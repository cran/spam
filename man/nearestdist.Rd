% This is file ../spam/man/nearestdist.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   



\name{nearestdist}
\alias{dist.spam}
\alias{nearest.dist}
\alias{distance}
\title{Distance Matrix Computation}
\description{This function computes and returns specific elements of distance
     matrix computed by
     using the specified distance measure.}
\usage{
nearest.dist( x, y=NULL, method = "euclidean",
             delta = 1, upper = if (is.null(y)) FALSE else NULL,
             p=2, miles=TRUE, R=NULL)
}
\arguments{
  \item{x}{Matrix of first set of locations where each row gives the
          coordinates of a particular point. See also \sQuote{Details}.}
  \item{y}{Matrix of second set of locations where each row gives the
          coordinates of a particular point. If this is missing \code{x} is
          used. See also \sQuote{Details}.}
  \item{method}{the distance measure to be used. This must be one of
          \code{"euclidean"}, \code{"maximum"}, \code{"minkowski"} or
                \code{"greatcircle"}. Any unambiguous substring can be
          given.}
  \item{delta}{only distances smaller than \code{delta} are recorded,
    see Details.} 
  \item{upper}{Should the entire matrix (\code{NULL}) or only the upper-triagonal (\code{TRUE})
               or lower-triagonal (\code{FALSE}) values be calculated.} 
  \item{p}{The power of the Minkowski distance.} 
  \item{miles}{For great circle distance: If true distances are in statute miles if false distances in
          kilometers.} 
  \item{R}{For great circle distance: Radius to use for sphere to find spherical distances. If \code{NULL}
          the radius is either in miles or kilometers depending on the
          values of the miles argument. If \code{R=1} then distances are of
          course in radians.}
%  \item{eps}{deprecated. Left for backwards consistency.} 
%  \item{diag}{deprecated. Left for backwards consistency. See \sQuote{Details}.} 
}
\value{A \code{spam} object containing the distances spanned between
  zero and \code{delta}. The sparse matrix may contain many zeros
  (e.g., collocated data). However, to calculate covariances, these zeros
  are essential.}

\details{For great circle distance, the matrices \code{x} and \code{y}
  contain the degrees longitudes in the first and the degrees latitudes
  in the second column.  \code{eps} and \code{delta} are in
  degrees. Hence to restrict to distances smaller than \code{delta.km},
  one has to specify \code{delta=delta.km*360/(6378.388*2*pi)}.

  The distance is
  in single precision (I am still not sure where I lose the double precision in
  the Fortran code) and if calculating the entire matrix
  \code{upper=NULL} (instead of adding its transpose) it may not
  pass the symmetry checks, for example.\cr 

  Default value of Earth's radius is 3963.34miles (6378.388km).\cr

  The formerly depreciated arguments \code{eps} and \code{diag} are now
  eliminated. 
  
  \code{x} and \code{y} can be any object with an existing
  \code{as.matrix} method.\cr

  A quick scan revealed distance functions in at least 7 packages
  (around 2008). The argument names should be as general as possible and
  be coherent with many (but not all) available distance functions.\cr
  
  The Fortran code is based on a idea of Doug Nychka.
    }
%\references{}
\seealso{\code{\link{spam_rdist}}}
\examples{
# Note that upper=T and using t(X)+X is quicker than upper=NULL;
#     upper=T marginally slower than upper=F.

# To compare nearest.dist with dist, use as.dist(...)
nx <- 4
x <- expand.grid(as.double(1:nx),as.double(1:nx))
sum( ( as.dist(nearest.dist( x, delta=nx*2))-
          dist(x)                            )^2)

# Create nearest neighbor structures:
par(mfcol=c(1,2))
x <- expand.grid(1:nx,1:(2*nx))
display( nearest.dist( x, delta=1))
x <- expand.grid(1:(2*nx),1:nx)
display( nearest.dist( x, delta=1))

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
