% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/display.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{display}
\alias{display}
\alias{display.spam}
\alias{display,spam-method}
%\alias{display,spam.chol.NgPeyton-method}
\title{Graphially Represent the Nonzero Entries}
\description{The function represents the nonzero entries in a simple
  bicolor plot.
}
\usage{
display(x, ...)
}
\arguments{\item{x}{matrix of class \code{spam} or \code{spam.chol.NgPeyton}.}
  \item{...}{any other arguments passed
    to \code{image.default}/\code{plot}.}
}
%\value{}
\details{\code{spam.getOption("imagesize")} determines  if the sparse matrix is
  coerced into a matrix and the plotted with \code{image.default} or if
  the matrix is simply represented as a scatterplot with
  \code{pch="."}.  The points are scaled according to
  \code{cex*getOption("spam.cex")/(nrow + ncol)}.
  For some devices or for non-square matrices,
  \code{cex} needs probably some adjustment.
}
%\references{}
\seealso{\code{\link{image}}, \code{\link{spam.options}}}
\examples{
set.seed(13)

smat <- spam_random(8)
par(mfcol=c(1,2), pty='s')
options(spam.imagesize = 1000)
display(smat)
options(spam.imagesize = 10)
display(smat, cex=.25)


# very large but very sparse matrix
smat <- spam_random(2^14, distribution=rnorm, density=1e-5, verbose=TRUE)
par(mfcol=c(1, 1), mai=c(.4,.4,.1,.1), pty='s')
display(smat)
}
\author{Reinhard Furrer}
\keyword{hplot}
