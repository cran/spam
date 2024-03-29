% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/triplet.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\name{triplet}
\alias{triplet}
\title{Transform a "spam" Format to Triplets}
\description{Returns a list containing the indices and elements of a  \code{spam} object.}
\usage{
triplet(x, tri=FALSE)
}
\arguments{
  \item{x}{sparse matrix of class \code{spam} or a matrix.}
  \item{tri}{Boolean indicating whether to create individual row and
    column  indices vectors.}
}
\value{A list with elements
  \item{indices}{a by two matrix containing the indices if
    \code{tri=FALSE}.}
  \item{i,j}{vectors  containing the row and column indices if \code{tri=TRUE}.}
  \item{values}{a vector containing the matrix elements.}
}
\details{
  The elements are row (column) first if \code{x} is a \code{spam} object (matrix).\cr
}
%\references{}
\seealso{ \code{\link{spam.creation}} for the inverse operation and
  \code{foreign} for other transformations.}
\examples{
x <- diag.spam(1:4)
x[2,3] <- 5
triplet(x)
all.equal( spam( triplet(x, tri=TRUE)), x)

}
\author{Reinhard Furrer}
\keyword{array}
