% This is file ../spam/man/bandwidth.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{bandwidth}
\alias{bandwidth}
\title{Bandwidth of a Sparse Matrix}
\description{Returns the lower and upper bandwidth of a sparse matrix}
\usage{
bandwidth(A)
}
\arguments{
  \item{A}{spam object}
}
\details{The matrix does not need to be diagonal. Values can be negative
  indicating the the matrix contains a band cinfined in the
  upper or lower triangular part.
	}
\value{Integer vector containing the lower and upper bandwidth}
%\references{}
\seealso{
	\code{\link{diag.spam}}.
}
\examples{
bandwidth(spam(c(0, 1), 3, 2))

bandwidth(spam(c(0, 0, 1, rep(0, 9)), 4, 3))

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
