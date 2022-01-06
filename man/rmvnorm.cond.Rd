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


\name{rmvnorm.conditional}
\alias{rmvnorm.cond}
\alias{rmvnorm.conditional}
\alias{rmvnorm.conditional.spam}
\title{Draw Conditional Multivariate Normals}
\description{
  Fast way to draw conditional multivariate normals when the covariance
  matrix is sparse.}
\usage{
rmvnorm.conditional(n, y, mu = rep.int(0, dim(SigmaXX)[1]+dim(SigmaYY)[1]),
                    SigmaXX, SigmaYY, SigmaXY, noise, RstructYY = NULL, ...)
}
\arguments{
\item{n}{number of observations.}
\item{y}{observed vector.}
\item{mu}{mean vector.}
\item{SigmaXX}{covariance of X, required (of class \code{spam}).}
\item{SigmaXY}{cross-covariance of X-Y, optional (of class \code{spam}).}
\item{SigmaYY}{covariance of Y, required (of class \code{spam}).}
\item{noise}{observational noice of Y, optional. See \sQuote{Details}. }
\item{RstructYY}{the Cholesky structure of \code{SigmaYY}.}
\item{\dots}{arguments passed to \code{chol}.}
}
\details{Quite often, we want to draw condional observations \deqn{X|y}
  from the model \deqn{Y=X+e}, where \deqn{X} has covariance matrix
  \code{SigmaXX} and \deqn{e} has white noise.\cr

  Covariance of \deqn{Y} can be specified by \code{SigmaYY} or
  \code{SigmaXX+diag(noise,)}. If \deqn{Y} and \deqn{X} do not have the
  same dimensions, \code{SigmaXY} needs to be specified. \cr


  
  The function also implmements a general multivariate model, where the
  we only observe part of the vector.  The components are first
  \deqn{X} then \deqn{Y}.
}
\seealso{\code{\link{rmvnorm.spam}}.
}
\examples{
set.seed(12)
N <- 300
y <- c(5, -5, -5, 5)
SigmaXX <- as.spam(.95^abs(outer(1:N,1:N,"-")), eps=1e-4)
sel <- c(10,100,120,300)        # where we observe y
SigmaXY <- SigmaXX[, sel]
SigmaYY <- SigmaXX[sel,sel] + diag.spam(.01, length(y)) # some noise
x <- rmvnorm.conditional(3, y, SigmaXX=SigmaXX, SigmaXY=SigmaXY,
                         SigmaYY=SigmaYY)
# unconditional sample:
ux <- rmvnorm(1, Sigma=SigmaXX)
matplot(t(rbind(x, ux)), type='l', lty=1)
points(sel, y, pch=19)  
}
\author{Reinhard Furrer}
\keyword{algebra}
