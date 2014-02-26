% This is file ../spam/man/circulant.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   



\name{circulant}
\alias{circulant.spam}
\title{Create Circulant Matrices}
\description{Creates a circulant matrix in \code{spam} format.}
\usage{
circulant.spam(x, n = NULL, eps = .Spam$eps)
}
\arguments{
  \item{x}{the first row to form the circulant matrix or a list
    containing the indices and the nonzero values.}
  \item{n}{if \code{x} is a list, the dimension of the matrix.}
  \item{eps}{A tolerance parameter:  elements of \code{x} such that \code{abs(x)
    <= eps} set to zero. Defaults to \code{eps = .Spam$eps}} 
}
\value{The circulant matrix in \code{spam} format.}
%\details{The vector \code{y} has to be of the same length as \code{x}
%  and its first element is discarded. 
%    }
%\references{}
\seealso{\code{\link[magic]{circulant}} from package \pkg{magic}, \code{\link{toeplitz.spam}}}
\examples{
circulant.spam(c(1,.25,0,0,0))
}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
