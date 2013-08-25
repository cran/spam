% This is file ../spam/man/adjacency.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{adjacency}
\alias{adjacency}
\alias{adjacency.spam}
\alias{adjacency.landkreis}
\alias{germany.graph}
\title{Administrative districts of Germany}
\description{Constructing the adjacency graph 
  of the administrative districts of Germany}
\usage{
adjacency.landkreis( loc)
}
\arguments{
\item{loc}{location of the graph structure, can be an URL.}
}
\details{The function is included as an example on how to
construct adjacency matrices form a (common) adjacency structure.
For the particular example, note that the nodes are not numbered
consecutively and that they start from zero.}
\value{a sparse matrix in \code{spam} format.}
\references{The adjacency data has been provided by Havard Rue and is
  also available in  \pkg{INLA}.}
\seealso{\code{\link{germany.plot}} super-seeding \code{map.landkreis}
  for plotting.\cr
  \code{\link{Oral}}.}
\examples{
\dontrun{
loc <- system.file("demodata/germany.adjacency", package="spam")
display( adjacency.landkreis( loc))
}
}
\author{Reinhard Furrer}
\keyword{hplot}
