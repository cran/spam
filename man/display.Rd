% This is file ../spam0.20-2/man/display.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{display}
\alias{display}
\alias{display.spam}
\alias{display,spam-method}
\alias{display,spam.chol.NgPeyton-method}
\title{Graphially represent the nonzero entries}
\description{The function represents the nonzero entries in a simple
  bicolor plot. 
}
\usage{
  display(x, ...)
}
\arguments{\item{x}{matrix of class \code{spam} or \code{spam.chol.NgPeyton}.}
%  \item{cex}{for very large matrices, the dot size may need to be scaled.}
  \item{...}{any other arguments passed
    to \code{image.default}/\code{plot}.}
}
%\value{}
\details{\code{spam.getOption('imagesize')} determines  if the sparse matrix is
  coerced into a matrix and the plotted with \code{image.default} or if
  the matrix is simply represented as a scatterplot with
  \code{pch="."}. The points are scaled according to
  \code{cex*spam.getOption('cex')/(nrow+ncol)}.
  For some devices or for non-square matrices,
  \code{cex} needs probably some adjustment.
}
%\references{}
\seealso{\code{\link{image}}, \code{\link{spam.options}}}
\examples{
set.seed(13)

nz <- 8
ln <- nz
smat <- spam(0,ln,ln)
smat[cbind(sample(ln,nz),sample(ln,nz))] <- 1:nz

par(mfcol=c(1,2),pty='s')
spam.options( imagesize=1000)
display(smat)
spam.options( imagesize=10)
display(smat)


# very large but very sparse matrix
nz <- 128
ln <- nz^2
smat <- spam(0,ln,ln)
smat[cbind(sample(ln,nz),sample(ln,nz))] <- 1:nz

par(mfcol=c(1,1),pty='s')
display(smat, cex=100)

}
\author{Reinhard Furrer}
\keyword{hplot}
