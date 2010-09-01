% This is file ../spam0.23-0/man/chol.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{chol}
\alias{chol.spam}
\alias{chol,ANY-method}
\alias{chol,matrix-method}
\alias{chol,spam-method}
\alias{update}
\alias{update,spam.chol.NgPeyton-method}
\alias{update.spam.chol.NgPeyton}
\alias{chol}              
\title{Cholesky Factorization for Sparse Matrices}
\description{
\code{chol} performs a Cholesky
decomposition of a symmetric positive definite sparse matrix \code{x} 
of class \code{spam}.}
\usage{
chol(x, \dots)
chol.spam(x, pivot = "MMD", method = 'NgPeyton', memory =
                 list(), eps =  .Spam$eps, \dots)

update.spam.chol.NgPeyton(object, x,...)
%chol(x, method, ordering, memory, \dots)
}
\arguments{
\item{x}{symmetric positive definite matrix of class \code{spam}.}
\item{pivot}{should the matrix be permuted, and if, with what
  algorithm, see Details below.}
\item{method}{Currently, only \code{NgPeyton} is implemented.}
\item{memory}{Parameters specific to the method, see Details below.}
\item{eps}{threshold to test symmetry. Defaults to \code{.Spam$eps}.}
\item{\dots}{further arguments passed to or from other methods.}
\item{object}{an object from a previous call to \code{chol}.}
}
\value{The function returns the Cholesky factor in an object of class
  \code{spam.chol.}\emph{method}. Recall that the latter is the Cholesky
  factor of a reordered matrix \code{x}, see also \code{\link{ordering}}. 
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

Often the sparsity structure is fixed and does not change, but the
entries do. In those cases, we can update the Cholesky factor with
\code{update.spam.chol.NgPeyton} by suppling a Cholesky factor and the
updated matrix.

The option \code{cholupdatesingular} determines how singular matrices
are handled by \code{update}. The function hands back an error
(\code{"error"}), a warning (\code{"warning"}) or the value \code{NULL}
(\code{"null"}).\cr



The Cholesky decompositions requires parameters, linked to memory
allocation.  If the default values are too small the Fortran routine
returns an error to \R, which allocates more space and calls the Fortran
routine again. The user can also pass better estimates of the allocation
sizes to \code{chol} with the argument \code{memory=list(nnzR=...,
  nnzcolindices=...)}. The minimal sizes for a fixed sparsity
structure can be obtained from a \code{summary} call.\cr


The output of \code{chol} can be used with \code{forwardsolve} and
\code{backsolve} to solve a system of linear equations.\cr


Notice that the Cholesky factorization of the package \code{SparseM} is also
based on the algorithm of Ng and Peyton (1993).  Whereas the Cholesky
routine of the package \code{Matrix} are based on 
\code{CHOLMOD} by Timothy A. Davis (\code{c} code).
}
\references{ 
Ng, E. G. and B. W. Peyton (1993) Block sparse Cholesky algorithms
on advanced uniprocessor computers, \emph{SIAM J. Sci. Comput.}, \bold{14},
1034--1056.

George, A. and  Liu, J. (1981)
\emph{Computer Solution of Large Sparse Positive Definite Systems},
Prentice Hall.
}
\note{
Although the symmetric structure of \code{x} is needed, only the upper
diagonal entries are used. By default, the code does check for
symmetry (contrarily to \code{base:::chol}). However,
depending on the matrix size, this is a time consuming test.
A test is ignored if
\code{.spam.options("cholsymmetrycheck")} is set to \code{FALSE}. 


If a permutation is supplied with \code{pivot},
\code{.spam.options("cholpivotcheck")} determines if the permutation is
tested for validity (defaults to \code{TRUE}). 
}
\seealso{\code{\link{det}}, \code{\link{solve}}, 
  \code{\link{forwardsolve}}, \code{\link{backsolve}} and \code{\link{ordering}}.
}
\examples{
# generate multivariate normals:
set.seed(13)
n <- 25    # dimension
N <- 1000  # sample size
Sigma <- .25^abs(outer(1:n,1:n,"-"))
Sigma <- as.spam( Sigma, eps=1e-4)

cholS <- chol( Sigma)    
# cholS is the upper triangular part of the permutated matrix Sigma
iord <- ordering(cholS, inv=TRUE)

R <- as.spam(cholS)
mvsample <- ( array(rnorm(N*n),c(N,n)) \%*\% R)[,iord]
# It is often better to order the sample than the matrix
# R itself. 

# 'mvsample' is of class 'spam'. We need to transform it to a
# regular matrix, as there is no method 'var' for 'spam' (should there?).
norm( var( as.matrix( mvsample)) - Sigma, type="HS")
norm( t(R) \%*\% R - Sigma, type="sup")


}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer, based on Ng and Peyton (1993) Fortran routines}
\keyword{algebra}
