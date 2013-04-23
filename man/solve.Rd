% This is file ../spam0.29-3/man/solve.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   







\name{spam solve}
\alias{solve.spam}
\alias{backsolve-methods}
\alias{backsolve,ANY-method}
\alias{backsolve,spam-method}
\alias{backsolve,matrix-method}
\alias{backsolve}
\alias{backsolve.spam}
\alias{forwardsolve-methods}
\alias{forwardsolve,ANY-method}
\alias{forwardsolve,spam-method}
\alias{forwardsolve,matrix-method}
\alias{forwardsolve}
\alias{forwardsolve.spam}
\alias{chol2inv}
\alias{chol2inv.spam}
\alias{chol2inv,spam-method}
\alias{chol2inv,spam.chol.NgPeyton-method}
\alias{solve,ANY-method}
\alias{solve,spam-method}
\title{Linear Equation Solving for Sparse Matrices}
\description{
  \code{backsolve} and \code{forwardsolve} solve a system
  of linear equations where the coefficient matrix
   is upper or lower triangular.
\cr
\code{solve} solves a linear system or computes the inverse
of a matrix if the right-hand-side is missing.
}
\usage{
# solve(a, b, ...)
# solve.spam(a, b, ...)

\S3method{solve}{spam}(a, b, ...)

# backsolve(r, x, ...)
backsolve.spam(r, x, ...)
# forwardsolve(l, x, ...)
forwardsolve.spam(l, x, ...)

chol2inv(x,...)
}
\arguments{
\item{a}{symmetric positive definite matrix of class \code{spam} or a Cholesky factor 
		  as the result of a \code{chol} call.}
\item{l,r}{object of class \code{spam} or \code{spam.chol.}\emph{method} returned by the function
\code{chol}.}
\item{x,b}{vector or regular matrix of right-hand-side(s) of a system of linear equations.}
\item{\dots}{further arguments passed to or from other methods, see \sQuote{Details} below.}
}
\details{
  We can solve \code{A \%*\% x = b} by first computing the Cholesky decomposition \code{A =
  t(R)\%*\%R)}, then solving \code{t(R)\%*\%y = b} for \code{y}, and
finally solving \code{R\%*\%x = y} for \code{x}.
   \code{solve} combines \code{chol}, a Cholesky decomposition of a
symmetric positive definite sparse matrix, with \code{forwardsolve} and
then \code{backsolve}.\cr

In case \code{a} is from a \code{chol}, then \code{solve} is an
   efficient way to calculate \code{backsolve(a, forwardsolve( t(a), b))}.


   However, for \code{a.spam} and \code{a.mat} from a \code{chol} call
   wit a sparse and ordinary matrix,  note that
   \code{forwardsolve( a.mat, b, transpose=T, upper.tri=T)}
   is equivalent to
   \code{forwardsolve( t(a.mat), b)}
   and \code{backsolve(a.spam, forwardsolve(a.spam, b, transpose=T, upper.tri=T))}
   yields the desired result. But
          \code{backsolve(a.spam,forwardsolve(t(a.spam), resid))} is
   wrong because  \code{t(a.spam)} is a \code{spam} and not a
   \code{spam.chol.NgPeyton} object.

   
\code{forwardsolve} and \code{backsolve} solve a system of linear
equations where the coefficient matrix is lower (\code{forwardsolve}) or
upper (\code{backsolve}) triangular.  Usually, the triangular matrix is
result from a \code{chol} call and it is not required to transpose it
for \code{forwardsolve}.  Note that arguments of the default
methods \code{k}, \code{upper.tri} and \code{transpose} do not have any
effects here.

Notice that it is more efficient to solve successively the linear
equations (both triangular solves) than to implement these in the
Fortran code. 

If the right-hand-side in \code{solve} is missing it will compute
the inverse of a matrix. For details about the specific Cholsesky
decomposition, see \code{\link{chol}}.

Recall that the Cholesky factors are from ordered matrices. 


\code{chol2inv(x)} is a faster way to \code{solve(x)}.
}
\note{There is intentionally no \acronym{S3} distinction between the classes
  \code{spam} and \code{spam.chol.}\emph{method}.}
\references{See references in \code{\link{chol}}. 
}
\seealso{\code{\link{det}}, \code{\link{chol}} and \code{\link{ordering}}.
}
\examples{
# Generate multivariate form a covariance inverse:
# (usefull for GRMF)
set.seed(13)
n <- 25    # dimension
N <- 1000  # sample size
Sigmainv <- .25^abs(outer(1:n,1:n,"-"))
Sigmainv <- as.spam( Sigmainv, eps=1e-4)


Sigma <- solve( Sigmainv)  # for verification 
iidsample <- array(rnorm(N*n),c(n,N))

mvsample <- backsolve( chol(Sigmainv), iidsample)
norm( var(t(mvsample)) - Sigma)

# compare with:
mvsample <- backsolve( chol(as.matrix( Sigmainv)), iidsample, n)
   #### ,n as patch 
norm( var(t(mvsample)) - Sigma)



# 'solve' step by step:
b <- rnorm( n)
R <- chol(Sigmainv)
norm( backsolve( R, forwardsolve( R, b))-
      solve( Sigmainv, b) )
norm( backsolve( R, forwardsolve( R, diag(n)))- Sigma )


}
% backsolve( chol(as.matrix(V)[ord,ord]),iidsample)[iord,]
%
\author{Reinhard Furrer, based on Ng and Peyton (1993) Fortran routines}
\keyword{algebra}
