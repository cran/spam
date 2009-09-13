% This is file ../spam0.15-5/man/complexity.Rd
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
\details{\code{chol} performs a Cholesky decomposition of a symmetric
positive definite sparse matrix \code{x} of class
\code{spam}. Currently, there is only the block sparse Cholesky
algorithm of Ng and Peyton (1993) implemented (\code{method=NgPeyton}).

To pivot/permute the matrix, you can choose between the multiple minimum
degree (\code{pivot=MMD}) or reverse Cuthill-Mckee (\code{pivot=RCM})
from George and Lui (1981). It is also possible to furnish a specific
permutation in which case \code{pivot} is a vector. For compatibility
reasons, \code{pivot} can also take a logical in which for \code{FALSE}
no permutation is done and for \code{TRUE} is equivalent to
\code{MMD}.\cr

Often the sparseness structure is fixed and does not change, but the
entries do. In those cases, we can update the Cholesky factor with
\code{update.spam.chol.NgPeyton} by suppling a Cholesky factor and the
updated matrix.

% link to here from a few instances.

The Cholesky decompositions requires parameters, linked to memory
allocation.  If the default values are too small the Fortran routine
returns an error to \R, which allocates more space and calls the Fortran
routine again. The user can also pass better estimates of the allocation
sizes to \code{chol} with the argument \code{memory=list(nnzR=...,
  nnzcolindices=...)}. The minimal sizes for a fixed sparseness
structure can be obtained from a \code{summary} call.\cr


The output of \code{chol} can be used with \code{forwardsolve} and
\code{backsolve} to solve a system of linear equations.\cr

}
\references{ 
Ng, E. G. and B. W. Peyton (1993) Block sparse Cholesky algorithms
on advanced uniprocessor computers, \emph{SIAM J. Sci. Comput.}, \bold{14},
1034--1056.

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
