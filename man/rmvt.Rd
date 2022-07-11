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


\name{rmvt}
\alias{rmvt}
\alias{rmvt.spam}
\title{Draw From a Multivariate t-Distribution}
\description{
  Fast ways to draw from a multivariate t-distribution the scale (covariance) matrix
  is sparse.}
\usage{
rmvt(n, Sigma, df = 1, delta = rep(0, nrow(Sigma)),
    type = c("shifted", "Kshirsagar"), ..., sigma)
rmvt.spam(n, Sigma, df = 1, delta = rep(0, nrow(Sigma)),
    type = c("shifted", "Kshirsagar"), ..., sigma)
}
\arguments{
\item{n}{number of observations.}
\item{Sigma}{scale matrix (of class \code{spam}).}
\item{df}{degrees of freedom.}
\item{delta}{vector of noncentrality parameters.}
\item{type}{type of the noncentral multivariate t distribution.}
\item{\dots}{arguments passed to \code{rmvnorm.spam}.}
\item{sigma}{similar to \code{Sigma}. Here for portability with \code{mvtnorm::rmvt()}}
}
\details{This function is very much like \code{rmvt()} from the package
  \pkg{mvtnorm}. We refer to the help of the afore mentioned.
}
\references{See references in \code{mvtnorm::rmvt()}.
}
\seealso{\code{\link{rmvnorm}}.
}
%\examples{}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer}
\keyword{algebra}
