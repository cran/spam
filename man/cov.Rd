% This is file ../spam/man/cov.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   

\name{covmat}
\alias{covmat}
\alias{cov.exp}
\alias{cov.sph}
\alias{cov.nug}
\alias{cov.wu1}
\alias{cov.wu2}
\alias{cov.wu3}
\alias{cov.wend1}
\alias{cov.wend2}
\alias{cov.mat}
\title{Covariance functions}
\description{
  Evaluate a covariance function.}
\usage{
covmat(h, theta, ... , type="sph")

cov.exp(h, theta, ... , eps= .Spam$eps)
cov.sph(h, theta, ... , eps= .Spam$eps)
cov.nug(h, theta, ... , eps= .Spam$eps)
cov.wu1(h, theta, ... , eps= .Spam$eps)
cov.wu2(h, theta, ... , eps= .Spam$eps)
cov.wu3(h, theta, ... , eps= .Spam$eps)
cov.wend1(h, theta, ... , eps= .Spam$eps)
cov.wend2(h, theta, ... , eps= .Spam$eps)
cov.mat(h, theta, ... , eps= .Spam$eps)

}
\arguments{
\item{h}{object containing the lags.}
\item{theta}{parameter of the covariance function, see
  \sQuote{Details}.}
\item{type}{covariance function specification.}
\item{\dots}{arguments passed from other methods.}
\item{eps}{tolerance level.}
}
\details{\code{covmat} is a wrapper that calls the other functions
  according to the argument \code{type}. The nomenclature is similar to
  \code{premat}
  \cr
  The parametrization is (range, sill, [smoothness], nugget), where
  only the range needs to be specified. Default values are (1,[1],0).
  In case of negative parameter values, a warning is issued and the
  absolute value is retained.
  Although more cryptic, having all arguments as a single vector
  simplifies optimization with \code{optim}. 
  \cr
  Currently, the functions distinguish between a sparse \code{spam}
  object \code{h} and any other numeric type. In the future, this might
  change and appropriate methods will be implemented. 
}
\value{Covariance function evaluated on \code{h}.}
%\note{There is intentionally no \acronym{S3} distinction between the classes
%  \code{spam} and \code{spam.chol.}\emph{method}.}
\references{Any classical book about geostatistics.}
\seealso{\code{\link{precmat}}.}
\examples{
locs <- cbind(runif(10),runif(10))
h <- nearest.dist(locs, delta=.3)
Sigma <- cov.sph(h, c(.3, 1, .1))

\dontrun{
h <- seq(0, to=1, length.out=100)
plot( h, cov.exp(h, c(1/3,1)), type='l', ylim=c(0,1))
type <- c("sph","wendland1","wendland2","wu1","wu2","wu3")
for (i in 1:6) 
  lines( h, covmat(h, 1, type=type[i]), col=i+1)
legend('topright',legend=type, col=2:7, lty=1)

}
}

\author{Reinhard Furrer}
\keyword{algebra}
