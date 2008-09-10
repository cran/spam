\name{allequal}
\alias{all.equal.spam}
\alias{all.equal,matrix,spam-method}
\alias{all.equal,spam,spam-method}
\title{Test if Two 'spam' Objects are (Nearly) Equal}
\description{Utility to compare two \code{spam} objects
     testing 'near equality'.  Depending on the type of difference, comparison is
     still made to some extent, and a report of the differences is
     returned.}
\usage{
all.equal.spam(target, current, tolerance = .Machine$double.eps^0.5,
    scale = NULL, check.attributes = FALSE,...)
}
\arguments{
  \item{target}{a \code{spam} object.}
  \item{current}{another \code{spam} object to be compared with \code{target}.}
  \item{tolerance}{numeric >= 0.  Differences smaller than
    \code{tolerance} are not considered.}
  \item{scale}{numeric scalar > 0 (or \code{NULL}). See 'Details'.}
  \item{check.attributes}{currently not yet implemented.}
  \item{...}{Further arguments for different methods.}
}
\value{Either \code{TRUE} or a vector of 'mode' \code{"character"} describing the
  differences between
  \code{target} and \code{current}.
}
\details{
  Numerical comparisons for
  \code{scale = NULL} (the default) are done by
     first computing the mean absolute difference of the two numerical
     vectors.  If this is smaller than \code{tolerance} or not finite,
     absolute differences are used, otherwise relative differences
     scaled by the mean absolute difference.

     If \code{scale} is positive, absolute comparisons are made after
     scaling (dividing) by \code{scale}.

     Don't use \code{all.equal.spam} directly in \code{if}
  expressions-either use \code{isTRUE(all.equal.spam(....))} or
  \code{identical} if appropriate.

  
  Cholesky decomposition routines use this function to test for
  symmetry.

  A method for \code{matrix-spam} objects is defined as well.
}
%\seealso{}
\examples{
obj <- diag.spam(2)
obj[1,2] <- .Machine$double.eps

all.equal( diag.spam(2), obj)

all.equal( t(obj), obj)

all.equal( t(obj), obj*1.1)

# We can compare a spam to a matrix
all.equal(diag(2),diag.spam(2))

# the opposite does often not make sense,
# hence, it is not implemented.
all.equal(diag.spam(2),diag(2))
}
\author{Reinhard Furrer}
\keyword{array}
