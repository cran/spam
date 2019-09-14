% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/landkreis.Rd.                       %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   










\name{map.landkreis}
\alias{map.landkreis}
\title{Administrative districts of Germany}
\description{Displaying data
  over the administrative districts of Germany}
\usage{
map.landkreis(data, col=NULL, zlim=range(data), add=FALSE,
              legendpos=c( 0.88,0.9,0.05,0.4))
}
\arguments{
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
\details{The function \code{\link{germany.plot}} super-seeds
  \code{map.landkreis} (it is several factors faster). \cr

The perfect position of the legend is an art per se and depends on
various \code{par} parameters. See also the source code of the function
\code{image.plot} of \pkg{fields}.}
\seealso{\code{\link{germany.plot}} super-seeding \code{map.landkreis}.}
\examples{
\dontrun{
data( Oral)
par( mfcol=c(1,2))
germany.plot( log( Oral$Y), legend=TRUE)
map.landkreis( log( Oral$Y))
}
}
\author{Reinhard Furrer}
\keyword{hplot}

% dev.off()
% dev.off();system.time( for (i in 1:20) map.landkreis(1:544))
% dev.off();system.time( for (i in 1:20) germany.plot(1:544))
