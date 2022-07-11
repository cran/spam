% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/print.Rd.                           %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{print}
\docType{methods}
\alias{print}
\alias{print.spam}
\alias{print,spam-method}
%\alias{print,spam.chol.NgPeyton-method}
\alias{print.spam.chol}
\alias{print.spam.chol.NgPeyton}
\alias{summary.spam}
\alias{summary,spam-method}
%\alias{summary,spam.chol.NgPeyton-method}
\alias{summary.spam.chol}
\alias{summary.spam.chol.NgPeyton}
\alias{print_nnzpos}
\title{Printing and Summarizing Sparse Matrices}
\description{Printing (non-zero elements) of sparse matrices
  and summarizing the sparsity structure thereof.
}
\usage{
  \S4method{print}{spam}(x, ...)
\S4method{summary}{spam}(object, ...)
}
\arguments{\item{x}{matrix of class \code{spam} or \code{spam.chol.}\emph{method}.}
  \item{object}{matrix of class \code{spam} or \code{spam.chol.}\emph{method}.}
  \item{...}{any other arguments passed to \code{print.default}.
             If the non-standard argument \code{minimal} is set to \code{FALSE},
             an extended spam print is available with logical argument \code{rowpointer}
             to print rowpointers, and \code{zerosymbol} defining the character
             to display the zero element.}
}
\value{\code{NULL} for \code{print}, because the information is printed with \code{cat} there
  is no real need to pass any object back. \cr
  %
  A list containing the non-zero elements and the density for
  \code{summary}
  for class \code{spam}.\cr
  %
  A list containing the non-zero elements of the factor, the density and
  the fill-in for  \code{summary} for class \code{spam.chol.NgPeyton}.}
\details{\code{getOption('spam.printsize')} determines  if the sparse matrix is
  coerced into a matrix and the printed as an array or if
  only the non-zero elements of the matrix are given.
}
%\references{}
\seealso{\code{\link{display}} or \code{\link{image}} for a graphical visualization;
   \code{\link{spam.options}}}
\examples{
set.seed(13)
smat <- spam_random(8)
par(mfcol=c(1,2),pty='s')
options(spam.printsize=1000)
print(smat)
options(spam.printsize=10)
print(smat)
summary(smat)
summary(smat)$nnz

smat@entries[1:5] <- 0
print(smat, minimal = FALSE)
print(smat, minimal = FALSE, rowpointer = TRUE)
smat@rowpointers
print_nnzpos(smat)

}
\author{Reinhard Furrer}
\keyword{hplot}
