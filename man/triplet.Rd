\name{triplet}
\alias{triplet}
\title{Transform a spam format to triplets}
\description{Returns a list containing the indices and elements of a  \code{spam} object.}
\usage{
triplet(x)
}
\arguments{
  \item{x}{sparse matrix of class \code{spam} or a matrix.}
}
\value{A list with elements 
  \item{indices}{a matrix containing the indices}
  \item{elements}{a vector containing the elements}
}
\details{
  The elements are row (column) first if \code{x} is a \code{spam} object (matrix).\cr
}
%\references{}
\seealso{ \code{\link{spam.list}} for the inverse operation and
  \code{foreign} for other transformations.}
\examples{
triplet(diag.spam(4))
}
\author{Reinhard Furrer}
\keyword{array}
