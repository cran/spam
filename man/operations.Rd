% This is file ../spam0.28-0/man/operations.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{spam operations}
\alias{Ops.spam}
\alias{spam.ops}
\alias{\%*\%-methods}
\alias{\%*\%,ANY,ANY-method}
\alias{\%*\%,spam,spam-method}
\alias{\%*\%,spam,matrix-method}
\alias{\%*\%,spam,numeric-method}
\alias{\%*\%,matrix,spam-method}
\alias{\%*\%,numeric,spam-method}
\alias{\%d*\%}
\alias{\%d*\%,spam,spam-method}
\alias{\%d*\%,spam,ANY-method}
\alias{\%d*\%,matrix,spam-method}
\alias{\%d*\%,matrix,ANY-method}
\alias{\%d*\%,spam,numeric-method}
\alias{\%d*\%,numeric,spam-method}
\alias{\%d*\%,numeric,matrix-method}
\alias{\%d*\%,numeric,numeric-method}
\alias{\%d+\%}
\alias{\%d+\%,spam,spam-method}
\alias{\%d+\%,spam,ANY-method}
\alias{\%d+\%,matrix,spam-method}
\alias{\%d+\%,matrix,ANY-method}
\alias{\%d+\%,spam,numeric-method}
\alias{\%d+\%,numeric,matrix-method}
\alias{\%d+\%,numeric,spam-method}
\alias{\%d+\%,numeric,numeric-method}
\alias{+,spam,spam-method}
\alias{+,ANY,spam-method}
\alias{+,spam,ANY-method}
\alias{-,spam,spam-method}
\alias{-,ANY,spam-method}
\alias{-,spam,ANY-method}
\alias{*,spam,spam-method}
\alias{*,ANY,spam-method}
\alias{*,spam,ANY-method}
\alias{/,spam,spam-method}
\alias{/,ANY,spam-method}
\alias{/,spam,ANY-method}
\alias{&,spam,ANY-method}
\alias{&,spam,spam-method}
\alias{&,ANY,spam-method}
\alias{|,spam,ANY-method}
\alias{|,ANY,spam-method}
\alias{|,spam,spam-method}
\alias{^,spam-method}
\alias{\%\%,spam-method}
\alias{\%/\%,spam-method}
\alias{>,spam-method}
\alias{>=,spam-method}
\alias{<,spam-method}
\alias{<=,spam-method}
\alias{==,spam-method}
\alias{!=,spam-method}
\alias{norm} 
\alias{norm.spam} 
\alias{norm,ANY-method} 
\alias{norm,spam,character-method} 
\alias{norm,spam,missing-method} 
%\alias{norm,numeric-method} 
%\alias{t,ANY-method} 
%\alias{t,spam-method} 
\alias{t.spam}
\alias{ncol,spam-method} 
\alias{nrow,spam-method} 
\alias{dim,ANY-method} 
\alias{dim,spam-method} 
\alias{[.spam}
\alias{[<-.spam}
\alias{[<-,spam,missing,missing-method}
\alias{[<-,spam,missing,vector-method}
\alias{[<-,spam,vector,missing-method}
\alias{[<-,spam,vector,vector-method}
\alias{[<-,spam,matrix,missing-method}
\alias{[<-,spam,matrix,matrix-method}
\alias{[<-,spam,spam,missing-method}
\alias{[<-,spam,ANY,ANY-method}
\alias{[,spam,ANY,ANY-method}
\alias{plot.spam}
\alias{subset.spam}
\alias{assign.spam}
\title{Basic Linear Algebra for Sparse Matrices}
\description{Basic linear algebra operations for sparse matrices
of class \code{spam}.  
}
%\usage{x \%*\% y
%y \%d*\% x
%y \%d+\% x
%x[i,]
%...
%}
%\arguments{
%\item{x}{matrix of class \code{spam}.}
%\item{y}{matrix of class \code{spam} or a dense matrix or vector.}
%\item{value}{replacement values.}
%\item{i,j}{vectors of elements to extract or replace.}
%\item{nrow}{optional number of rows for the result.}
%}
\details{Linear algebra operations for matrices of class 
\code{spam} are designed to behave exactly as for 
regular matrices.  In particular, matrix multiplication, transpose, addition, 
subtraction and various logical operations should work as with the conventional
dense form of matrix storage, as does indexing, rbind, cbind, and diagonal
assignment and extraction (see for example \code{\link{diag}}).
Further functions with identical behavior are \code{dim} and thus \code{nrow}, \code{ncol}.




The function \code{norm} calculates the (matrix-)norm of the argument.
The argument \code{type} specifies the \code{l1} norm, \code{sup} or max
norm (default), or the Frobenius or Hilbert-Schmidt
(\code{frobenius/hs}) norm. Partial matching can be used. For example,
\code{norm} is used to check for symmetry in the function \code{chol} by
computing the norm of the difference between the matrix and its
transpose

The operator \code{\%d*\%} efficiently multiplies a diagonal matrix (in
vector form) and a sparse matrix and is used for compatibility with the
package fields. More specifically, this method is used in the internal
functions of \code{Krig} to make the code more readable. It avoids
having a branch in the source code to handle the diagonal or nondiagonal
cases.  Note that this operator is not symmetric: a vector in
the left argument is interpreted as a diagonal matrix and a vector in
the right argument is kept as a column vector.

The operator \code{\%d+\%} efficiently adds a diagonal matrix (in vector
form) and a sparse matrix, similarly to the operator \code{\%d+\%}.

}
\references{Some Fortran functions are based on 
  \url{http://www-users.cs.umn.edu/~saad/software/SPARSKIT/sparskit.html}
}
\seealso{
\code{\link{spam}} for coercion and other class relations involving the
sparse matrix classes.
}
\examples{
# create a weight matrix and scale it:
\dontrun{
wij <- distmat
# with distmat from a nearest.dist(..., upper=TRUE) call

n <- dim(wij)[1]

wij@entries <- kernel( wij@entries, h) # for some function kernel
wij <- wij + t(wij) + diag.spam(n)     # adjust from diag=FALSE, upper=TRUE

sumwij <- wij \%*\% rep(1,n)
    # row scaling:
    #   wij@entries <- wij@entries/sumwij[ wij@colindices]
    # col scaling:
wij@entries <- wij@entries/sumwij[ rep(1:n, diff(wij@rowpointers))]
}
}
\keyword{algebra}
