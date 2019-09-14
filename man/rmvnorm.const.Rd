% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/rmvnorm.const.Rd.                   %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{rmvnorm.const}
\alias{rmvnorm.const}
\alias{rmvnorm.prec.const}
\alias{rmvnorm.canonical.const}
\title{Draw Constrainted Multivariate Normals}
\description{
  Fast ways to draw multivariate normals with linear constrains when the variance or precision matrix
  is sparse.}
\usage{
rmvnorm.const(n, mu = rep.int(0, dim(Sigma)[1]), Sigma, Rstruct = NULL,
              A = array(1, c(1,dim(Sigma)[1])), a=0, U=NULL,  ...)
rmvnorm.prec.const(n, mu = rep.int(0, dim(Q)[1]), Q, Rstruct = NULL,
              A = array(1, c(1,dim(Q)[1])), a=0, U=NULL,  ...)
rmvnorm.canonical.const(n, b, Q, Rstruct = NULL,
              A = array(1, c(1,dim(Q)[1])), a=0, U=NULL, ...)
}
\arguments{
\item{n}{number of observations.}
\item{mu}{mean vector.}
\item{Sigma}{covariance matrix of class \code{spam}.}
\item{Q}{precision matrix.}
\item{b}{vector determining the mean.}
\item{Rstruct}{the Cholesky structure of \code{Sigma} or \code{Q}.}
\item{A}{Constrain matrix.}
\item{a}{Constrain vector.}
\item{U}{see below.}
\item{\dots}{arguments passed to \code{chol}.}
}
\details{The functions \code{rmvnorm.prec} and \code{rmvnorm.canonical}
  do not requrie sparse precision matrices.
  For \code{rmvnorm.spam}, the differences between regular and sparse
  covariance matrices are too significant to be implemented here.  \cr
  Often (e.g., in a Gibbs sampler setting), the sparsity structure of
  the covariance/precision does not change. In such setting, the
  Cholesky factor can be passed via \code{Rstruct} in which only updates
  are performed (i.e., \code{update.spam.chol.NgPeyton} instead of a
  full \code{chol}).
}
%\note{There is intentionally no \acronym{S3} distinction between the classes
%  \code{spam} and \code{spam.chol.}\emph{method}.}
\references{See references in \code{\link{chol}}.
}
\seealso{\code{\link{rmvnorm.spam}}.
}
\examples{
# to be filled in

}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer}
\keyword{algebra}
