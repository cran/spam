% This is file ../spam0.29-3/man/germany.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{germany}
\alias{map.landkreis}
\alias{adjacency.landkreis}
\alias{germany.graph}
\title{Administrative districts of Germany}
\description{Constructing the adjacency graph and displaying the data
  over the administrative districts of Germany}
\usage{
adjacency.landkreis( loc)
map.landkreis(data, col=NULL, zlim=range(data), add=FALSE,
              legendpos=c( 0.88,0.9,0.05,0.4))
}
\arguments{
\item{loc}{location of the graph structure, can be an URL.}
\item{data}{vector of length 544}
\item{col}{color scheme to be used. By default uses \code{tim.colors} if
available or a generic gray scale.}
\item{zlim}{the minimum and maximum values for which colors should be
          plotted, defaulting to the range of \code{data}.}
\item{add}{logical, if true adds to current plot.}
\item{legendpos}{if package \pkg{fields} is loaded, puts a legend at that position.}
}
\references{The code of \code{map.landkreis} is very similar to
\code{germany.map} from the package \pkg{INLA}.}
\details{\code{adjacency.landkreis} is included as an example on how to
construct adjacency matrices form a (common) adjacency structure.
For the particular example, note that the nodes are not numbered
consecutively and that they start from zero.\cr
The perfect position of the legend is an art per se and depends on
various \code{par} parameters. See also the source code of the function
\code{image.plot} of \pkg{fields}.}
\seealso{\code{\link{Oral}}.}
\examples{
loc <- system.file("demodata/germany.adjacency", package="spam")
display( adjacency.landkreis( loc))

map.landkreis( Oral$E)
}
\author{Reinhard Furrer}
\keyword{hplot}
