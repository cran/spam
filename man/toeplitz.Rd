% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/toeplitz.Rd.                        %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{toeplitz}
\alias{toeplitz.spam}
\title{Create Toeplitz Matrices}
\description{Creates symmetric and asymmetric Toeplitz matrices.}
\usage{
toeplitz.spam(x, y = NULL, eps = getOption("spam.eps"))
}
\arguments{
  \item{x}{the first row to form the Toeplitz matrix.}
  \item{y}{for asymmetric Toeplitz matrices, this contains the first
    column.}
  \item{eps}{A tolerance parameter:  elements of \code{x} such that \code{abs(x)
    <= eps} set to zero. Defaults to \code{eps = getOption("spam.eps")}.} 
}
\value{The Toeplitz matrix in \code{spam} format.}
\details{The vector \code{y} has to be of the same length as \code{x}
  and its first element is discarded. 
    }
%\references{}
\seealso{\code{\link{toeplitz}}, \code{\link{circulant.spam}}}
\examples{
toeplitz.spam(c(1,.25,0,0,0))

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
