% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/cov.Rd.                             %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{covmat}
\alias{covmat}
\alias{cov.exp}
\alias{cov.sph}
\alias{cor.sph}
\alias{cov.nug}
\alias{cov.wu1}
\alias{cov.wu2}
\alias{cov.wu3}
\alias{cov.wend1}
\alias{cov.wend2}
\alias{cov.mat}
\alias{cov.mat12}
\alias{cov.mat32}
\alias{cov.mat52}
\alias{cov.finnmat}
\title{Covariance Functions}
\description{
  Evaluate a covariance function.}
\usage{
covmat(h, theta, ... , type="sph")

cov.exp(h, theta, ... , eps= getOption("spam.eps"))
cov.sph(h, theta, ... , eps= getOption("spam.eps"))
cov.nug(h, theta, ... , eps= getOption("spam.eps"))
cov.wend1(h, theta, ... , eps= getOption("spam.eps"))
cov.wend2(h, theta, ... , eps= getOption("spam.eps"))
cov.wu1(h, theta, ... , eps= getOption("spam.eps"))
cov.wu2(h, theta, ... , eps= getOption("spam.eps"))
cov.wu3(h, theta, ... , eps= getOption("spam.eps"))
cov.mat(h, theta, ... , eps= getOption("spam.eps"))
cov.finnmat(h, theta, ... , eps= getOption("spam.eps"))
cov.mat12(h, theta, ... , eps= getOption("spam.eps"))
cov.mat32(h, theta, ... , eps= getOption("spam.eps"))
cov.mat52(h, theta, ... , eps= getOption("spam.eps"))

cor.sph(h, range, ... , eps= getOption("spam.eps"))

}
\arguments{
\item{h}{object containing the lags.}
\item{theta}{parameter of the covariance function, see
  \sQuote{Details}.}
\item{range}{parameter defining the compact support.}
\item{type}{covariance function specification.}
\item{\dots}{arguments passed from other methods.}
\item{eps}{tolerance level, see
  \sQuote{Details}.}
}
\details{\code{covmat} is a wrapper that calls the other functions
  according to the argument \code{type}. The nomenclature is similar to
  \code{precmat}.
  \cr
  The parametrization is (range, [partial-sill = 1], [smoothness = 1], [nugget = 0]), where
  only the range needs to be specified.
  In case of negative parameter values, a warning is issued and the
  absolute value is retained.
  Although more cryptic, having all arguments as a single vector
  simplifies optimization with \code{optim}.
  \cr
  The parameters are and locations are up to precision \code{epsilon}.
  That means that all distances smaller than \code{eps} are considered zero; a nugget smaller than \code{eps} is ignored; a range smaller than \code{eps} represents a nugget model; etc.
  \cr
  \code{cov.finnmat()} is similar to \code{cov.mat()} but with the
  \code{sqrt(8*smoothness)/range} argument in the Bessel function
  (instead of \code{1/range}).
  \code{cov.mat12()} is a wrapper to \code{cov.exp()}
  \code{cov.mat32()}, and \code{cov.mat52()} are fast version of
  \code{cov.mat()} with smoothness 3/2 and 5/2, respectively (factor 10).
  \cr
  \code{cor.sph(,range)} is a fast version of \code{cov.sph(,c(range,1,0))}.
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
set.seed(123)
n <- 200
locs <- cbind(runif(n),runif(n))
h <- nearest.dist(locs, delta=sqrt(2), upper = NULL)
Sigma <- cov.sph(h, c(.3, 1, .1))

iidsample <- rnorm(n)
cholS <- chol.spam(as.spam(Sigma))
iorder <- iord <- ordering(cholS, inv = TRUE)
sample <- (iidsample \%*\% as.spam(cholS))[iorder]
plot(locs, col = fields::tim.colors(n = 256)[cut(sample, n)], pch = 20)

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
