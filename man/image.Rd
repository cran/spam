% This is file ../spam0.29-0/man/image.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{image}
\alias{image}
\alias{image.spam}
\alias{image,spam-method}
%\alias{image,spam.chol.NgPeyton-method}
\title{Display a spam Object as Color Image}
\description{The function creates a grid of colored rectangles with
  colors corresponding to the values in \code{z}. }
\usage{
image(x, ...)
}
\arguments{\item{x}{matrix of class \code{spam} or \code{spam.chol.NgPeyton}.}
%  \item{cex}{for very large matrices, the dot size may need to be scaled.}
  \item{...}{any other arguments passed
    to \code{image.default} and \code{plot}.}
}
%\value{}
\details{\code{getOption("imagesize")} determines  if the sparse matrix is
  coerced into a matrix and the plotted similarly to \code{image.default} or if
  the matrix is simply represented as a scatterplot with
  \code{pch="."}. The points are scaled according to
  \code{cex*spam.getOption("cex")/(nrow+ncol)}.
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
nz <- 8
ln <- nz
smat <- spam(0, ln, ln)
smat[ cbind(sample(ln,nz), sample(ln,nz))] <- 1:nz

par(mfcol=c(1,2),pty='s')
spam.options( imagesize=1000)
image(smat) # better: col=tim.colors(nz))
spam.options( imagesize=10)
image(smat) # better: col=tim.colors(nz))


nz <- 128
ln <- nz^2
smat <- spam(0,ln,ln)
smat[cbind(sample(ln,nz), sample(ln,nz))] <- 1:nz

par(mfcol=c(1,1), pty='s')
image(smat,cex=100) # better:, col=tim.colors(nz))

}
\author{Reinhard Furrer}
\keyword{hplot}
