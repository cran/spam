% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/germanydata.Rd.                     %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   





\name{germany.data}
\alias{germany.data}
\alias{germany.info}
\alias{germany.poly}
\alias{germany.spam}
\alias{germany}
\docType{data}

\title{Meta-data about administrative districts of Germany}
\description{Data for the display of  data
  over the administrative districts of Germany}
\format{\code{germany.info} is a list with elements 
  \describe{
    \item{n}{544 (number of districts around 1990).}
    \item{xrep,yrep}{representative coordinates of the districts (vectors
      of length 544)}
    \item{xlim,ylim}{2-vectors defining the limits of the districts.}
    \item{polyid}{linking the polygons to the districts (599 vector).}
    \item{id}{linking the districts to Community Identification Number.}}
  \code{germany.poly} defines the polygons. It is a 17965 by two matrix,
  each polygon separated by a row of \code{NA}s, each district by two rows.\cr
  \code{germany} defines the polygons in form of a list (backwards compatibility). }


\references{The meta-data has been constructed based on (essentially)
  files from the   package \pkg{INLA}, see \code{demo(Bym)}.\cr

  See also \url{http://de.wikipedia.org/wiki/Amtlicher_Gemeindeschl\%C3\%BCssel}
and \url{http://de.wikipedia.org/wiki/Liste_der_Landkreise_in_Deutschland}}
\details{The representative coordinates are calculated based on the mean
  value of the polygon coordinates. This creates sometimes strange
  values, e.g., district Leer.}
\seealso{\code{\link{germany.plot}} \code{\link{Oral}}.}
\examples{
# Plot the Bundeslaender:
germany.plot(germany.info$id\%/\%1000,col=rep(2:8,3), legend=FALSE)
}
\author{Reinhard Furrer}
\keyword{hplot}
