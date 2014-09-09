% This is file ../spam/man/ordering.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   

\name{ordering}
\docType{methods}
\alias{ordering}
\alias{ordering.spam}
\alias{ordering-methods}
%\alias{ordering,spam.chol.NgPeyton-method}
\alias{ordering.spam.chol}
\alias{ordering.spam.chol.NgPeyton}
\alias{ordering,spam-method}
\alias{ordering,matrix-method}
\title{Extract the permutation}
\description{Extract the (inverse) permutation used by the Cholesky
  decomposition}
\usage{
ordering( x, inv=FALSE)
}
\arguments{
  
\item{x}{object of class \code{spam.chol.}\emph{method} returned by the function
\code{chol}.}
\item{inv}{Return the permutation (default) or inverse thereof.}
}
\details{
  Recall that calculating a Cholesky factor from a sparse matrix
  consists of finding a permutation first, then calculating the factors
  of the permuted matrix. The ordering is important when working with
  the factors themselves.\cr

  The ordering from a full/regular matrix is \code{1:n}.\cr
  
  Note that there exists many different algorithms to find
  orderings.  
  \cr
  
  See the examples, they speak more than 10 lines.
}
\seealso{\code{\link{chol.spam}}, \code{\link{solve.spam}}.
}
\examples{
# Construct a pd matrix S to work with (size n)
n <- 100    # dimension
S <- .25^abs(outer(1:n,1:n,"-"))
S <- as.spam( S, eps=1e-4)
I <- diag(n)  # Identity matrix

cholS <- chol( S)
ord <- ordering(cholS)
iord <- ordering(cholS, inv=TRUE)

R <- as.spam( cholS ) # R'R = P S P', with P=I[ord,],
  # a permutation matrix (rows permuted).
RtR <- t(R) \%*\% R

# the following are equivalent:
as.spam( RtR -            S[ord,ord] )
as.spam( RtR[iord,iord] - S )
as.spam( t(R[,iord]) \%*\% R[,iord] - S )

# trivially:
as.spam( t(I[iord,]) - I[ord,])  # (P^-1)' = P  
as.spam( t(I[ord,]) - I[,ord])  # 
as.spam( I[iord,] - I[,ord])
as.spam( I[ord,]\%*\%S\%*\%I[,ord] - S[ord,ord] )
   # pre and post multiplication with P and P' is ordering
}
\author{Reinhard Furrer}
\keyword{algebra}
