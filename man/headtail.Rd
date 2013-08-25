% This is file ../spam/man/headtail.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{head}
\alias{head}
\alias{head.spam}
\alias{head,spam-method}
\alias{head,spam.chol.NgPeyton-method}
\alias{tail}
\alias{tail.spam}
\alias{tail,spam-method}
\alias{tail,spam.chol.NgPeyton-method}
\title{
  Return the First or Last Part of an Object
}
\description{
  Returns the upper left or lower right part of a
  \code{\linkS4class{spam}} object.
}
\usage{
head(x, \dots)
tail(x, \dots)
# head.spam(x, n = 6L, m = n, \dots)
# tail.spam(x, n = 6L, m = n, addrownums = TRUE, \dots)
}
\arguments{
  \item{x}{a \code{\linkS4class{spam}} object}
  \item{n}{a single integer. If positive, size for the resulting
    object: number of elements for a vector (including lists), rows for
    a matrix or data frame or lines for a function. If negative, all but
    the \code{n} last/first number of elements of \code{x}.}
  \item{m}{similar to \code{n} but for the number of columns.}
  \item{addrownums}{create row and column namves them from the selected
    elements.}
  \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
  For matrices, 2-dim tables and data frames, \code{head()} (\code{tail()}) returns
  the first (last) \code{n} rows and \code{m} columns when \code{n > 0} or all but the
  last (first) \code{n} rows when \code{n < 0} (with similar behavior
  for \code{m}).

  \code{tail()} will add row and column names of
  the form \code{"[n,]"} and \code{"[,n]"} to the result, so that it looks similar to the
  last lines and columns of \code{x} when printed.  Setting \code{addrownums =
    FALSE} suppresses this behaviour.


  A method for \code{\linkS4class{spam.chol.NgPeyton}} objects is
  exported as well.
}
\value{
  An regular matrix.
}
\author{
  Reinhard Furrer
}
\examples{
head( precmat.RW2( 10))
tail( precmat.season(n=10, season=3), n=4, m=10)
}
\keyword{ manip }
