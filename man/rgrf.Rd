% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/rmvnorm.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{rgrf}
\alias{rgrf}
\alias{rgrf.spam}
\title{Draw From a Gaussian Random Field}
\description{
  Fast and intuitive ways to draw from a Gaussian random field.}
\usage{
rgrf( n,
   locs, nx, ny=nx, xlim=c(0,1), ylim=c(0,1), tau=0,
   Covariance, theta, beta=0, X, 
   method=c('chol'),  method.args=list(sparse=FALSE), 
   eps = getOption("spam.eps"), drop=TRUE, attributes=TRUE, ...)
}
\arguments{
\item{n}{number of observations.}
\item{locs}{locations, the result of \code{as.matrix(locs)} will be
  used.}
\item{nx,ny}{if no locations are specified, at least one of these to
  specify the grid dimension.}
\item{xlim,ylim}{Domain, see \sQuote{Details}.}
\item{tau}{perturbation degree, see \sQuote{Details}.}
\item{Covariance}{covariance function name.}
\item{theta}{covariance parameter.}
\item{beta}{mean or vector for regression-type mean.}
\item{X}{design matrix for regression-type mean.}
 \item{method}{based on Choleski factorization.}
 \item{method.args}{list of arguments that can be passed to the corresponding approach.
  For \code{"chol"} it can be, e.g., \code{RStruct}, \code{chol.args}, \code{cov.args}.}
\item{eps}{small value, anything smaller is considered a collocation.}
\item{drop}{logical, if a single realization should be returned as a vector.}
\item{attributes}{logical, if should attributes be passed back.}

\item{\dots}{currently not used.}
}
\details{If no locations are given, the function constructs these
  according a regular or a regular perturbed grid. The perturbation is
  determined by \code{tau}, which has to be greater than zero (no
  perturbation) and strictly smaller than 1/2 (max perturbation).

  The regular grid has spacing (here for x) \code{dx=diff(xlim)/nx} and runs
  from \code{xlim[1]+dx/2} to \code{xlim[2]-dx/2}.
  The locations are at least  \code{(1/nx-2*tau*dx)} separated.
  \cr
  
  Currently, the only method implemented is a Cholesky factorization
  routine,  (much as in \code{\link{rmvnorm}}).

  The \code{rdist()} from the \code{fields} package is awefully
  fast. Unless one has very sparse covariance matrices, a sparse
  approach is not bringing a lot of improvements. 

  
  The methods may use different covariance construction approaches and
  thus the nesting of \code{cov.args} in \code{method.args}.
}
\seealso{\code{\link{rgrf}}, \code{\link{chol}} and \code{\link{ordering}}.
}
\examples{
require(fields)
# Regular grid with constant mean:
nx <- 10
field <- rgrf(1, nx=nx,  Covariance="cov.wend2", theta=c(.5, 1), beta=5)
quilt.plot(cbind(attr(field,"locs"),z=field), nx=nx, ny=nx)
points(attr(field,"locs"))

# Irregluar grid:
field <- rgrf(1, nx=10, tau=0.3, Covariance="cov.mat", theta=c(.2, 1, 1.5))
fields::quilt.plot(attr(field,"locs"), field)

}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer}
\keyword{algebra}
