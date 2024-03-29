% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/USprecip.Rd.                        %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{USprecip}
\alias{USprecip}
\docType{data}
\title{Monthly Total Precipitation (mm) for April 1948 in the Contiguous United States}
\description{
This is a useful spatial data set of moderate to large size consisting of 11918
locations. See \url{https://www.image.ucar.edu/GSP/Data/US.monthly.met/} for the source of these data.
}

\format{
This data set is an array containing the following columns:

\describe{
\item{lon,lat}{
  Longitude-latitude position of monitoring stations.
}
\item{raw}{
  Monthly total precipitation in millimeters for April 1948.
}
\item{anomaly}{
  Preipitation anomaly for April 1948.
}
\item{infill}{
Indicator, which station values were observed (5906 out of the 11918)
compared to which were estimated.
}
}
}
\source{\url{https://www.image.ucar.edu/GSP/Data/US.monthly.met/}}
\seealso{\code{\link[fields]{RMprecip}}}
\references{
Johns, C., Nychka, D., Kittel, T., and Daly, C. (2003)
Infilling sparse records of spatial fields.
\emph{Journal of the American Statistical Association},
98, 796--806.
}
\examples{
# plot
\dontrun{
library(fields)

data(USprecip)
par(mfcol=c(2,1))
quilt.plot(USprecip[,1:2],USprecip[,3])
US( add=TRUE, col=2, lty=2)
quilt.plot(USprecip[,1:2],USprecip[,4])
US( add=TRUE, col=2, lty=2)
}
}
\keyword{datasets}
