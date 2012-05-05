% This is file ../spam0.29-1/man/complexity.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{complexity}
\alias{complexity}
\alias{complexities}      
\title{Complexity for Sparse Matrices}
\description{A few results of computational complexities for selected
  sparse algoritms in \code{spam}
}
\details{A Cholesky factorization of an n-matrix requires n^3/3 flops. In case of
  banded matrices (bandwidth p, p<<n) a factorization requires about 2np^2
  flops. Forward- and backsolves for banded matrices require essentially
  2np flops.

  George and Liu (1981) proves that any reordering would require at
  least O(n^3/2) flops for the factorization and produce at least O(n
  log(n)) fill-ins for square lattices with a local neighbor hood.\cr
  They also show that algorithms based on nested dissection are optimal
  in the order of magnitude sense.

  More to follow.

  



%  It is also possible to furnish a specific
%permutation in which case \code{pivot} is a vector. For compatibility
%reasons, \code{pivot} can also take a logical in which for \code{FALSE}
%no permutation is done and for \code{TRUE} is equivalent to
%\code{MMD}.\cr%

%The exact operation count is often very difficult to determine. Easier
%bounds are obtained in the case of banded matrices.

}
\references{ 
%Ng, E. G. and B. W. Peyton (1993) Block sparse Cholesky algorithms
%on advanced uniprocessor computers, \emph{SIAM J. Sci. Comput.}, \bold{14},
%1034--1056.
%
George, A. and  Liu, J. (1981)
\emph{Computer Solution of Large Sparse Positive Definite Systems},
Prentice Hall.
}
\seealso{\code{\link{det}}, \code{\link{solve}}, 
  \code{\link{forwardsolve}}, \code{\link{backsolve}} and \code{\link{ordering}}.
}
%\examples{
%}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\keyword{algebra}
