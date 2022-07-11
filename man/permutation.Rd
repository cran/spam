% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/permutation.Rd.                     %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


\name{permutation}
\alias{permutation}
\alias{permutation.spam}
\alias{permutation,spam-method}
\alias{permutation,matrix-method}
\title{Permute a Matrix}
\description{Row and/or column permutes a (sparse) matrix.
}
\usage{
permutation.spam(A, P=NULL, Q=NULL, ind=FALSE, check=TRUE)
}
\arguments{
\item{A}{sparse matrix}
\item{P}{vector giving the row permutation.}
\item{Q}{vector giving the column permutation.}
\item{ind}{are the indices given. See examples.}
\item{check}{Should rudimentary checks be performed.}
}
\value{A permuted matrix.
}
\details{If P and Q are permutation matrices, the result is
  PAQ. However, it is also possible to specify the indices
  and to perform in a very efficient way
  \code{A[rowind, colind]}, see examples.


  A row permutation is much faster than a colum permutation.
  For very large matrices, a double transpose might be faster.

  The spam option \code{spam.checkpivot} determines if the permutation
  is verified.
}
%\references{
%}
\seealso{\code{\link{ordering}}, \code{\link{spam.options}}.}
\examples{
A <- spam(1:12,3)
P <- c(3,1,2)
Q <- c(2,3,1,4)

permutation(A,P,Q)-A[order(P),order(Q)]

permutation(A,P,Q,ind=TRUE)-A[P,Q]

}
\author{Reinhard Furrer}
\keyword{array}