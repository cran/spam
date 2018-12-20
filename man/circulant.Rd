% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/circulant.Rd.                       %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      %
%    Roman Flury [ctb], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{circulant}
\alias{circulant.spam}
\title{Create Circulant Matrices}
\description{Creates a circulant matrix in \code{spam} format.}
\usage{
circulant.spam(x, n = NULL, eps = getOption("spam.eps"))
}
\arguments{
  \item{x}{the first row to form the circulant matrix or a list
    containing the indices and the nonzero values.}
  \item{n}{if \code{x} is a list, the dimension of the matrix.}
  \item{eps}{A tolerance parameter:  elements of \code{x} such that \code{abs(x)
    <= eps} set to zero. Defaults to \code{eps = getOption("spam.eps")}} 
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
