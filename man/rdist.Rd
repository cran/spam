% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/rdist.Rd.                           %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      %
%    Roman Flury [ctb], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{fields wrapper}
\alias{spam_rdist}
\alias{rdist.dist}
\alias{spam_rdist.earth}
\title{Wrapper for Distance Matrix Computation}
\description{These functions are simple wrappers to \code{nearest.dist}
  to be used in \code{fields}.}
\usage{
spam_rdist( x1, x2, delta = 1)
spam_rdist.earth( x1, x2, delta = 1, miles=TRUE, R=NULL)
}
\arguments{
  \item{x1}{Matrix of first set of locations where each row gives the
          coordinates of a particular point. }
  \item{x2}{Matrix of second set of locations where each row gives the
          coordinates of a particular point.}
  \item{delta}{only distances smaller than \code{delta} are recorded,
    see Details.} 
  \item{miles}{For great circle distance: If true distances are in statute miles if false distances in
          kilometers.} 
  \item{R}{Radius to use for sphere to find spherical distances. If \code{NULL}
          the radius is either in miles or kilometers depending on the
          values of the miles argument. If \code{R=1} then distances are of
          course in radians.}
}
\value{A \code{spam} object containing the distances spanned between
  zero and \code{delta}. The sparse matrix may contain many zeros
  (e.g., collocated data). However, to calculate covariances, these zeros
  are essential.}
\details{These functions are wrappers to \code{rdist} and
  \code{rdist.earth} in \code{fields}. They are used
  to simplify the use of sparse matrices in functions like \code{mKrig}.
  \cr

  For great circle distance, the matrices \code{x1} and \code{x2}
  contain the degrees longitudes in the first and the degrees latitudes
  in the second column.  \code{delta} is in degrees. Hence to restrict
  to distances smaller than \code{delta.km}, one has to specify
  \code{delta=delta.km*360/(6378.388*2*pi)}.

}
%\references{}
\seealso{\code{\link{nearest.dist}}}
\examples{
\dontrun{
require(fields)
look <- mKrig(x,Y, Covariance="Wendland", dimension=2, k=1,
    cov.args=list( Distance='spam_rdist'))
}

}
\author{Reinhard Furrer}
\keyword{array}
