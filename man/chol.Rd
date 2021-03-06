% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/chol.Rd.                            %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{chol}
\alias{chol}
\alias{chol.spam}
\alias{chol,ANY-method}
\alias{chol,matrix-method}
\alias{chol,spam-method}
\alias{update.spam}
\alias{update,spam.chol.NgPeyton-method}
\alias{update.spam.chol.NgPeyton}
\title{Cholesky Factorization for Sparse Matrices}
\description{
\code{chol} performs a Cholesky
decomposition of a symmetric positive definite sparse matrix \code{x}
of class \code{spam}.}
\usage{
# chol(x, \dots)

\S4method{chol}{spam}(x, pivot = "MMD", method = "NgPeyton",
       memory = list(), eps = getOption("spam.eps"), Rstruct=NULL,
       \dots, verbose=FALSE)
%force64 =  getOption("spam.force64"),\dots)

# update.spam.chol.NgPeyton(object, x,...)
\S4method{update}{spam.chol.NgPeyton}(object, x,...)
%chol(x, method, ordering, memory, \dots)
}
\arguments{
\item{x}{symmetric positive definite matrix of class \code{spam}.}
\item{pivot}{should the matrix be permuted, and if, with what
  algorithm, see \sQuote{Details} below.}
\item{method}{Currently, only \code{NgPeyton} is implemented.}
\item{memory}{Parameters specific to the method, see \sQuote{Details} below.}
\item{eps}{threshold to test symmetry. Defaults to \code{getOption("spam.eps")}.}
%\item{force64}{logical vector of length 1. If \code{TRUE}, a 64-bit
%    spam matrix is returned in any case. If \code{FALSE}, a 32-bit
%    matrix is returned when possible. }
\item{Rstruct}{sparsity structure  of the factor, see \sQuote{Details} below.}
\item{\dots}{further arguments passed to or from other methods.}
\item{object}{an object from a previous call to \code{chol}, i.e.,
  sparsity structure of the factor.}
\item{verbose}{provides more details about the decomposition. Useful
  when working with  huge matrices.}
}
\value{The function returns the Cholesky factor in an object of class
  \code{spam.chol.}\emph{method}. Recall that the latter is the Cholesky
  factor of a reordered matrix \code{x}, see also \code{\link{ordering}}.
}

\details{\code{chol} performs a Cholesky decomposition of a symmetric
positive definite sparse matrix \code{x} of class
\code{spam}. Currently, there is only the block sparse Cholesky
algorithm of Ng and Peyton (1993) implemented (\code{method="NgPeyton"}).

To pivot/permute the matrix, you can choose between the multiple minimum
degree (\code{pivot="MMD"}) or reverse Cuthill-Mckee (\code{pivot="RCM"})
from George and Lui (1981). It is also possible to furnish a specific
permutation in which case \code{pivot} is a vector. For compatibility
reasons, \code{pivot} can also take a logical in which for \code{FALSE}
no permutation is done and for \code{TRUE} is equivalent to
\code{MMD}.\cr

Often the sparsity structure is fixed and does not change, but the
entries do. In those cases, we can update the Cholesky factor with
\code{update.spam.chol.NgPeyton} by suppling a Cholesky factor and the
updated matrix. For \code{U <- chol(A)},  \code{update(U,
  Anew)} and \code{chol(Anew, Rstruct=U)} are equivalent.

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
structure can be obtained from a \code{summary} call, see \sQuote{Examples}.\cr


The output of \code{chol} can be used with \code{forwardsolve} and
\code{backsolve} to solve a system of linear equations.\cr


Notice that the Cholesky factorization of the package \code{SparseM} is also
based on the algorithm of Ng and Peyton (1993).  Whereas the Cholesky
routine of the package \code{Matrix} are based on
\code{CHOLMOD} by Timothy A. Davis (\code{C} code).
}
\references{
Ng, E. G. and Peyton, B. W.  (1993) Block sparse Cholesky algorithms on
advanced uniprocessor computers, \emph{SIAM J. Sci. Comput.}, \bold{14},
1034--1056.


Gilbert, J. R., Ng, E. G. and Peyton, B. W. (1994) An efficient
   algorithm to compute row and column counts for sparse Cholesky
   factorization, \emph{SIAM J. Matrix Anal. Appl.}, \bold{15},
   1075--1091.


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
\code{options("spam.cholsymmetrycheck")} is set to \code{FALSE}.


If a permutation is supplied with \code{pivot},
\code{options("spam.cholpivotcheck")} determines if the permutation is
tested for validity (defaults to \code{TRUE}).
}
\seealso{\code{\link{det.spam}}, \code{\link{solve.spam}},
  \code{\link{forwardsolve.spam}}, \code{\link{backsolve.spam}} and \code{\link{ordering}}.
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
norm( var( as.matrix( mvsample)) - Sigma, type='m')
norm( t(R) \%*\% R - Sigma)


# To speed up factorizations, memory allocations can be optimized:
opt <- summary(cholS)
# here, some elements of Sigma may be changed...
cholS <- chol( Sigma, memory=list(nnzR=opt$nnzR,nnzcolindices=opt$nnzc))

}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer, based on Ng and Peyton (1993) Fortran routines}
\keyword{algebra}
