% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/kronecker.Rd.                       %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{kronecker}
\alias{kronecker.default}
\alias{kronecker.spam}
\alias{kronecker,spam,ANY-method}
\alias{kronecker,ANY,spam-method}
\alias{kronecker,spam,spam-method}
\title{Kronecker Products on Sparse Matrices}
\description{Computes the generalised kronecker product of two arrays, \code{X} and
     \code{Y}.}
\usage{
kronecker.spam(X, Y, FUN = "*", make.dimnames = FALSE, ...)
}
\arguments{
  \item{X}{sparse matrix of class \code{spam}, a vector or a matrix.}
  \item{Y}{sparse matrix of class \code{spam}, a vector or a matrix.}
  \item{FUN}{a function; it may be a quoted string. See details}
  \item{make.dimnames}{Provide dimnames that are the product of the dimnames of
    \code{X} and \code{Y}. }
    \item{ ...}{optional arguments to be passed to \code{FUN}.}
}
\value{An array \code{A} with dimensions
     \code{dim(X) * dim(Y)}.}
\details{The sparsity structure is determined by the ordinary
  \code{\%x\%}. Hence, the result of \code{kronecker(X, Y, FUN = "+")} is
    different depending on the input.
    }
%\references{}
%\seealso{\code{\link{chol}}}
\examples{
# Starting with non-spam objects, we get a spam matrix
kronecker.spam( diag(2), array(1:4, c(2, 2)))

kronecker( diag.spam(2), array(1:4, c(2, 2)))

# Notice the preservation of sparsity structure:
kronecker( diag.spam(2), array(1:4, c(2, 2)), FUN="+")


}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
