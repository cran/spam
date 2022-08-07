% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/spam.creation.Rd.                   %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{spam}
\alias{spam.creation}
\alias{initialize,spam-method}
\alias{spam}
\alias{spam.list}
\alias{spam.numeric}
\alias{spam,list-method}
\alias{spam,numeric-method}
%\alias{spam.spam}
%\alias{spam,spam-method}
\alias{as.spam,matrix-method}
\alias{as.spam,numeric-method}
\alias{as.spam,spam-method}
\alias{as.spam,dist-method}
\alias{as.spam,list-method}
\alias{as.spam}
\alias{as.spam.spam}
\alias{as.spam.numeric}
\alias{as.spam.matrix}
\alias{as.spam.chol.NgPeyton}
\alias{as.spam.dist}
\alias{as.spam.list}
\alias{is.spam}
\title{Sparse Matrix Class}
\description{
This group of functions evaluates and coerces changes in class structure.
}
\usage{
spam(x, nrow = 1, ncol = 1, eps = getOption("spam.eps"))
%     force64 = getOption("spam.force64"))

as.spam(x, eps = getOption("spam.eps"))
%, force64 = getOption("spam.force64"))

is.spam(x)
}
\value{
  A valid \code{spam} object.\cr
  \code{is.spam} returns \code{TRUE} if \code{x} is a \code{spam} object.}
\arguments{
  \item{x}{is a matrix (of either dense or sparse form), a list, vector
    object or a distance object}
\item{nrow}{number of rows of matrix }
\item{ncol}{number of columns of matrix }
%\item{force64}{logical vector of length 1. If \code{TRUE}, a 64-bit
%  spam matrix is returned in any case. If \code{FALSE}, a 32-bit
%  matrix is returned when possible. }
\item{eps}{A tolerance parameter:  elements of \code{x} such that
  \code{abs(x) < eps} set to zero. Defaults to \code{eps =
    getOption("spam.eps")} }
}
\details{
  The functions \code{spam} and \code{as.spam} act like \code{matrix}
  and \code{as.matrix}
  to coerce an object to a sparse matrix object of class \code{spam}.


  If \code{x} is a list, it should contain either two or three elements.
  In case of the former, the list should contain a \code{n} by two
  matrix of indicies (called \code{ind}) and the values.
  In case of the latter, the list should contain three vectors
  containing the row, column indices (called \code{i} and
  \code{j}) and the values. In both cases partial matching is done.
  In case there are several triplets with the same \code{i}, \code{j},
  the values are added.

  \code{eps} should be at least as large as \code{.Machine$double.eps}.


  If \code{getOption("spam.force64")} is \code{TRUE}, a 64-bit
  spam matrix is returned in any case. If \code{FALSE}, a 32-bit
  matrix is returned when possible.
}
\note{The zero matrix has the element zero stored in (1,1).\cr

  The functions do not test the presence of \code{NA/NaN/Inf}. Virtually
  all call a Fortran routine with the \code{NAOK=NAOK}
  argument, which defaults to \code{FALSE} resulting in an error.
  Hence, the \code{NaN} do not always properly propagate through (i.e.
  \code{spam} is not IEEE-754 compliant).

}
\references{Reinhard Furrer, Stephan R. Sain (2010).
  "spam: A Sparse Matrix R Package with Emphasis on MCMC
  Methods for Gaussian Markov Random Fields.",
        \emph{Journal of Statistical Software}, 36(10), 1-25,
        \doi{10.18637/jss.v036.i10}.
}
\seealso{
  \code{\link{SPAM}} for a general overview of the package;
  \code{\link{spam_random}} to create matrices with a random sparsity pattern;
  \code{\link{cleanup}} to purge a sparse matrix;
  \code{\link{spam.options}} for details about the \code{safemode} flag;
    \code{\link{read.MM}} and \code{\link{foreign}} to create \code{spam}
	matrices from MatrixMarket
    files and from certain \pkg{Matrix} or \pkg{SparseM} formats.
}
\examples{
# old message, do not loop, when you create a large sparse matrix
set.seed(13)
nz <- 128
ln <- nz^2
smat <- spam(0,ln,ln)
is <- sample(ln,nz)
js <- sample(ln,nz)
## IGNORE_RDIFF_BEGIN
system.time(for (i in 1:nz) smat[is[i], js[i]] <- i)
system.time(smat[cbind(is,js)] <- 1:nz)
## IGNORE_RDIFF_END

getClass("spam")


options(spam.NAOK=TRUE)
as.spam(c(1, NA))
}
\author{Reinhard Furrer}
\keyword{algebra}
