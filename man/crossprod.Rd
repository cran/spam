% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/crossprod.Rd.                       %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{crossprod}
\alias{crossprod.default}
\alias{crossprod.spam}
\alias{tcrossprod.spam}
\alias{crossprod,spam,missing-method}
\alias{tcrossprod,spam,missing-method}
\alias{crossprod,ANY,spam-method}
\alias{tcrossprod,ANY,spam-method}
\alias{crossprod,spam,spam-method}
\alias{tcrossprod,spam,spam-method}
\alias{crossprod,spam,ANY-method}
\alias{tcrossprod,spam,ANY-method}
\title{Spam Matrix Crossproduct}
\description{
     Given matrices \code{x} and \code{y} as arguments, return a matrix
     cross-product.  This is formally equivalent to (but usually
     slightly faster than) the call \code{t(x) \%*\% y} (\code{crossprod.spam}) or \code{x
     \%*\% t(y)} (\code{tcrossprod.spam}).}

\usage{

     crossprod.spam(x, y = NULL, \dots)
     
     tcrossprod.spam(x, y = NULL, \dots)
}
\arguments{
  \item{x, y}{matrices: \code{y = NULL} is taken to be the
          same matrix as \code{x}.  Vectors are promoted to single-column or
          single-row matrices, depending on the context.}
  \item{\dots}{potentially further arguments from other methods.}
}
\value{A double matrix}
\note{
     When \code{x} or \code{y} are not matrices, they are treated as column or
     row matrices.
   }
%\references{}
%\seealso{\code{\link{chol}}}
\examples{
crossprod.spam(diag.spam(2),1:2)

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
