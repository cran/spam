% This is file ../spam0.28-0/man/apply.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{apply}
\alias{apply.spam}
\title{Apply Functions Over Sparse Matrix Margins}
\description{Returns a vector or array or list of values obtained by applying a
     function to margins of a sparse matrix.}
\usage{apply.spam(X, MARGIN=NULL, FUN, ...)
}
\arguments{
  \item{X}{the \code{spam} matrix to be used.}
  \item{MARGIN}{a vector giving the subscripts which the function will be
          applied over. \code{1} indicates rows, \code{2} indicates columns,
          \code{NULL} or \code{c(1,2)} indicates rows and columns.}
  \item{FUN}{the function to be applied.}
  \item{...}{optional arguments to \code{FUN}.}
}

\details{This is a handy wrapper to apply a function to the (nonzero)
  elements of a sparse matrix. 
  For example, it is possible to apply a covariance matrix to a distance
  matrix obtained by \code{nearest.dist}, see Examples.\cr

  A call to \code{apply} only coerces the sparse matrix to a regular one.\cr
  
  The basic principle is applying the function to \code{@entries}, or to
  the extracted columns or rows (\code{[,i,drop=F]} or
  \code{[i,,drop=F]}). It is important to note that an empty column
  contains at least one zero value and may lead to non intuitive
  results.\cr

  This function may evolve over the next few releases. 
}
\value{Similar as a call to \code{apply} with a regular matrix. The most
  important cases are as follows. The
  result is a vector (\code{MARGIN} is length 1 and \code{FUN} is
  scalar) or a matrix (\code{MARGIN} is length 1 and \code{FUN} returns
  fixed length vectors, or \code{MARGIN} is length 2 and \code{FUN} is
  scalar) or a list (if \code{FUN} returns vectors of different lengths).}
%\references{}
\seealso{\code{base:apply} for more details on Value.}
\examples{
S <- as.spam(dist(1:5))
S <- apply.spam(S/2, NULL, exp)
# instead of 
# S@entries <- exp( S@entries/2) 

# Technical detail, a null matrix consists
# of one zero element.
apply.spam(S,c(1,2),pmax)
apply.spam(S,1,range)


# A similar example as for the base apply.
# However, no dimnames else we would get warnings. 
x <- as.spam(cbind(x1 = 3, x2 = c(0,0,0, 5:2)))
apply.spam(x, 2, mean, trim = .2)
col.sums <- apply.spam(x, 2, sum)
row.sums <- apply.spam(x, 1, sum)
rbind(cbind(x, row.sums), c(col.sums, sum(col.sums)))

apply.spam(x, 2, is.vector)

# Sort the columns of a matrix
# Notice that the result is a list due to the different
# lengths induced by the nonzero elements
apply.spam(x, 2, sort)

# Function with extra args:
cave <- function(x, c1, c2) c(mean(x[c1]), mean(x[c2]))
apply(x,1, cave,  c1=1, c2=c(1,2))

ma <- spam(c(1:4, 0, 0,0, 6), nrow = 2)
ma
apply.spam(ma, 1, table)  #--> a list of length 2
apply.spam(ma, 1, stats::quantile)# 5 x n matrix with rownames
}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
