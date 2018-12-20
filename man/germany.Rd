% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/germany.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      %
%    Roman Flury [ctb], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   





\name{germany.plot}
\alias{germany.plot}

\title{Plot administrative districts of Germany}
\description{Displaying data
  over the administrative districts of Germany}
\usage{
germany.plot(vect,  col=NULL, zlim=range(vect), legend=TRUE, 
             main=NULL, cex.axis=1, cex.main=1.5, add=FALSE, ... )
}
\arguments{
\item{vect}{vector of length 544}
\item{col}{color scheme to be used. By default uses \code{colorRampPalette(brewer.pal(9,"Blues"))(100)}.}
\item{zlim}{the minimum and maximum values for which colors should be
  plotted, defaulting to the range of \code{data}.}
\item{legend}{Should the legend be added, see \sQuote{Details}.}
\item{main}{an overall title for the plot.}
\item{cex.axis}{label size of legend.}
\item{cex.main}{label size of overall plot title.} 
\item{add}{logical, if true adds to current plot.}
\item{\dots}{additional arguments passed to \code{polygon}.}
}
\references{See also \url{http://de.wikipedia.org/wiki/Amtlicher_Gemeindeschl\%C3\%BCssel}
and \url{http://de.wikipedia.org/wiki/Liste_der_Landkreise_in_Deutschland}
%The code of \code{map.landkreis} is very similar to
%\code{germany.map} from the package \pkg{INLA}.}
}
\details{The legend is only added, provided (a) function
  \code{image.plot} is available.\cr

  The perfect position of the legend is an art per se and depends on
  various \code{par} parameters. One possiblity for finer control is to
  not plot it and to manually call the function
\code{image.plot} of \pkg{fields}.}
\seealso{\code{\link{Oral}}.}
\examples{
data( Oral)
germany.plot( Oral$Y/Oral$E)


# Plot the Bundeslaender:
germany.plot(germany.info$id\%/\%1000,col=rep(2:8,3), legend=FALSE)
}
\author{Reinhard Furrer}
\keyword{hplot}
