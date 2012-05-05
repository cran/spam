% This is file ../spam0.29-1/man/rowcolstats.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{rowSums}
\alias{rowSums.spam}
\alias{colSums.spam}
\alias{rowMeans.spam}
\alias{colMeans.spam}
\alias{rowSums}
\alias{colSums}
\alias{rowMeans}
\alias{colMeans}
\alias{rowSums,spam-method}
\alias{colSums,spam-method}
\alias{rowMeans,spam-method}
\alias{colMeans,spam-method}
\title{
Form Row and Column Sums and Means
}
\description{
Form row and column sums and means for sparse \code{\linkS4class{spam}} matrices
}
\usage{
rowSums(x, na.rm = FALSE, dims = 1, \dots)
colSums(x, na.rm = FALSE, dims = 1, \dots)
rowMeans(x, na.rm = FALSE, dims = 1, \dots)
colMeans(x, na.rm = FALSE, dims = 1, \dots)
}
\arguments{
  \item{x}{a \code{\linkS4class{spam}} object}
  \item{na.rm}{currently ignored}
  \item{dims}{ignored as we have only two dimensions.}
  \item{\dots}{potentially further arguments from other methods.}
}
\details{
Depending on the flag \code{}.
}
\value{
Vector of appropriate length.
}
%\references{}
\author{
Reinhard Furrer
}
%\note{}
\seealso{
\code{\link{apply.spam}}, \code{\link{spam.options}}.
}
\examples{
x <- spam( rnorm(20), 5, 4)
rowSums( x)
c( x \%*\% rep(1,4))
}
\keyword{manip}
