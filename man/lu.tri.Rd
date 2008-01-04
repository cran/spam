\name{triagonal}
\alias{lower.tri}
\alias{upper.tri}
\alias{lower.tri.spam}
\alias{upper.tri.spam}
\title{Lower and Upper Triangular Part of a Sparse Matrix}
\description{Returns the lower or upper triangular structure or entries of  a sparse matrix.
}
\usage{
lower.tri(x, diag = FALSE)
upper.tri(x, diag = FALSE)
}
\arguments{
  \item{x}{a sparse matrix of class \code{spam}}
  \item{diag}{logical. Should the diagonal be included?}
}
\details{Often not only the structure of the matrix is required but the
  entries as well. For compatibility, the default is only a structure
  consisting of ones (representing \code{TRUE}s). Setting the flag
  \code{spam.getOption('trivalues')} to \code{TRUE},
%  \code{.Spam$trivalues} to \code{TRUE},
  the function returns the actual values.
  }
\seealso{\code{\link{spam.options}} and \code{\link{diag}}}
\examples{
smat <- spam( c( 1,2,0,3,0,0,0,4,5),3)
upper.tri( smat)
upper.tri( smat, diag=TRUE)

%.Spam$trivalues <- TRUE
spam.options(trivalues=TRUE)
upper.tri( smat)

}
\keyword{array}
\keyword{algebra}
