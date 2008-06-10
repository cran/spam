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
             eps = .Spam$eps, delta = 1,
             diag = FALSE, upper = NULL,
             p=2, miles=TRUE, R=NULL)
}
\arguments{
  \item{x}{Matrix of first set of locations where each row gives the
          coordinates of a particular point. See also Details.}
  \item{y}{Matrix of second set of locations where each row gives the
          coordinates of a particular point. If this is missing \code{x} is
          used. See also Details.}
  \item{method}{the distance measure to be used. This must be one of
          \code{"euclidean"}, \code{"maximum"}, \code{"minkowski"} or
                \code{"greatcircle"}. Any unambiguous substring can be
          given.}
  \item{eps}{distances smaller than this number are considered zero.} 
  \item{delta}{only distances smaller than \code{delta} are recorded.} 
  \item{diag}{Should the diagonal be included? Only of importance if \code{y} is missing.} 
  \item{upper}{Should the entire matrix (\code{NULL}) or only the upper-triagonal (\code{TRUE})
               or lower-triagonal (\code{FALSE}) values be calculated.} 
  \item{p}{The power of the Minkowski distance.} 
  \item{miles}{If true distances are in statute miles if false distances in
          kilometers.} 
  \item{R}{Radius to use for sphere to find spherical distances. If \code{NULL}
          the radius is either in miles or kilometers depending on the
          values of the miles argument. If \code{R=1} then distances are of
          course in radians.}
}
\value{A \code{spam} object containing the distances spanned by \code{eps} and \code{delta}.}

\details{For the great circle distance, \code{x} and \code{y} contain
  the longitudes in the first and the latitudes in the second column.
  \code{x} and \code{y} can be any object with an existing
  \code{as.matrix} method.\cr

  A quick scan revieled distance functions in at least 7 packages. The
  argument names should be as general as possible and be coherend with
  many (but not all) available distance functions.\cr
  The Fortran code is based on a idea of Doug Nychka.
    }
%\references{}
%\seealso{\code{\link{chol}}}
\examples{
# Note that upper=T and using t(X)+X is quicker than upper=NULL;
#     upper=T marginally slower than upper=F.

# To compare nearest.dist with dist, use diag=FALSE, upper=TRUE
nx <- 4
x <- expand.grid(as.double(1:nx),as.double(1:nx))
sum( (nearest.dist( x, delta=nx*2, diag=FALSE, upper=TRUE)@entries-
              c(dist(x)))^2)

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
