\name{Summary}
\alias{Summary}
\alias{Summary,spam-method}
\alias{max}
\alias{min}
\alias{range}
\alias{prod}
\alias{sum}
\alias{any}
\alias{all}
\title{Summary functions}
\description{Applies the summary functions to 'spam' objects
}
\usage{max(x, ..., na.rm = FALSE)
min(x, ..., na.rm = FALSE)
range(x, ..., na.rm = FALSE)
prod(x, ..., na.rm = FALSE)
sum(x, ..., na.rm = FALSE)
any(x, ..., na.rm = FALSE )
all(x, ..., na.rm = FALSE)
}
\arguments{\item{x}{spam object.}
\item{...}{numeric arguments.}
\item{na.rm}{a logical indicating whether missing values should be
          removed, see details.}
}
\value{All functions operate on the vector \code{x@entries} and return the
  result thereof.
}
\details{The for this generic class typical \code{na.rm} argument has
    no weight here as NA/NaN/Inf are not meaningful (yet).\cr


%  Is implemented for R>=2.3.x only. However,
%  it would be possible to use Martin's proposed workaround:
%    \url{http://tolstoy.newcastle.edu.au/R/help/05/12/18192.html}
}
%\references{
%}
%\seealso{
%}
\examples{
smat <- diag.spam(rnorm(15))
range(smat)
}
\author{Reinhard Furrer}
\keyword{manip}