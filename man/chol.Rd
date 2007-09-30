\name{chol}
\alias{chol.spam}
\alias{chol,ANY-method}
\alias{chol,matrix-method}
\alias{chol,spam-method}              
\alias{chol}              
\alias{chol.default}              
\title{Cholesky Factorization for Sparse Matrices}
\description{
\code{chol} performs a Cholesky
decomposition of a symmetric positive definite sparse matrix \code{x} 
of class \code{spam}.}
\usage{
chol(x, \dots)
%chol(x, method, ordering, memory, \dots)
}
\arguments{
\item{x}{symmetric positive definite matrix of class \code{spam}.}
%\item{method}{Currently, only \code{NgPeyton} is implemented.}
%\item{ordering}{Currently, only \code{MinDeg} is implemented.}
%\item{memory}{Parameters specific to the method, see Details below.}
\item{\dots}{further arguments passed to or from other methods.}
}
\value{The function returns the Cholesky factor in an object of class
  \code{spam.chol.}\emph{method}. Recall that the latter is the Cholesky
  factor of a reordered matrix \code{x}, see also \code{\link{ordering}}. 
}

  \details{\code{chol} performs a Cholesky decomposition of a symmetric
positive definite sparse matrix \code{x} of class
\code{spam}. Currently, there is only the block sparse Cholesky
algorithm of Ng and Peyton (1993) implemented (\code{method=NgPeyton}),
based on a minimum degree permutation (\code{ordering=MinDeg}).\cr The
Cholesky routine is essentially identitcal to what is available in the
package \code{SparseM}.\cr

The output of \code{chol} can be used with 
\code{forwardsolve} \code{backsolve} to solve a system of linear equations.\cr

The Cholesky decompositions may require additional parameters, linked to
memory allocation. 
These memory  parameters and their default values for the  \code{NgPeyton}
method are \code{tmpmax=500*nrow(x)}, \code{nsubmax=length(x)},
\code{nnzlmax=max(4*nsubmax, floor(.2*nsubmax^1.3))}. In practice, only
the first one needs to be set to a larger value if one calculates the
Cholesky factor from a virtually full matrix.


}
\references{ 
Ng, E. G. and B. W. Peyton (1993), "Block sparse Cholesky algorithms
on advanced uniprocessor computers", \emph{SIAM J. Sci. Comput.}, \bold{14},
pp. 1034-1056.

Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}
}
\seealso{\code{\link{det}}, \code{\link{solve}}, 
  \code{\link{forwardsolve}}, \code{\link{backsolve}} and \code{\link{ordering}}.
}
\examples{
# generate multivariate normals:
n <- 25    # dimension
N <- 1000  # sample size
Sigma <- .25^abs(outer(1:n,1:n,"-"))
Sigma <- as.spam( Sigma, eps=1e-4)

cholS <- chol( Sigma)    
# cholS is the upper triangular part of the permutated matrix Sigma
iord <- ordering(cholS, inv=TRUE)

R <- as.spam(cholS)[,iord]
mvsample <- array(rnorm(N*n),c(N,n)) \%*\% R

# 'mvsample' is of class 'spam'. We need to transform it to a
# regular matrix, as there is no method 'var' for 'spam' (should there?).
norm( var( as.matrix( mvsample)) - Sigma, type="HS")
norm( t(R) \%*\% R - Sigma, type="sup")


}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer, based on Ng and Peyton (1993) Fortran routines}
\keyword{algebra}
