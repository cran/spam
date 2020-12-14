% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/image.Rd.                           %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{image}
\alias{image}
\alias{image.spam}
\alias{image.spam.chol}
\alias{image,spam-method}
%\alias{image,spam.chol.NgPeyton-method}
\title{Display a spam Object as Color Image}
\description{The function creates a grid of colored rectangles with
  colors corresponding to the values in \code{z}. }
\usage{
\S4method{image}{spam}(x, cex = NULL, ...)
}
\arguments{\item{x}{matrix of class \code{spam} or \code{spam.chol.NgPeyton}.}
  \item{cex}{for very large matrices, the dot size may need to be scaled.}
  \item{...}{any other arguments passed
    to \code{image.default} and \code{plot}.}
}
%\value{}
\details{\code{getOption("spam.imagesize")} determines  if the sparse matrix is
  coerced into a matrix and the plotted similarly to \code{image.default} or if
  the matrix is simply represented as a scatterplot with
  \code{pch="."}. The points are scaled according to
  \code{cex*getOption("spam.cex")/(nrow+ncol)}.
  For some devices or for non-square matrices,
  \code{cex} needs probably some adjustment.\cr
  The a zero matrix in \code{spam} format has as (1,1) entry the value
  zero and only missing entries are interpreted as \code{NA} in the
  scatter plot.
}
%\references{}
\seealso{\code{\link{display}} and  \code{\link{spam.options}}.\cr
  The code is based on  \code{\link[graphics]{image}} of \code{graphics}.
}
\examples{
set.seed(13)

smat <- spam_random(8)

par(mfcol=c(1,2),pty='s')
options(spam.imagesize=1000)
image(smat) # or use better color schemes
options(spam.imagesize=10)
image(smat, cex=.25)

smat <- spam_random(2^14, distribution=rnorm, density=1e-5, verbose=TRUE)
par(mfcol=c(1,1), mai=c(.4,.4,.1,.1), pty='s')
image(smat)

}
\author{Reinhard Furrer}
\keyword{hplot}
