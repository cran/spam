\name{.Spam}
\alias{.Spam}
\alias{spam-parameters}
\title{spam specific parameters}
\description{\code{.Spam} is a list holding information used for
  constructing and displaying spam objects.

}
\usage{.Spam
}
%\arguments{}
\value{\item{eps=.Machine$double.eps}{values smaller than this are
  considered as zero. This is only used when creating spam objects.}

\item{drop=FALSE}{default parameter for \code{drop} when subsetting}

\item{printsize=100}{the max number of elements of a matrix which we
  display as regular matrix.}

\item{imagesize=10000}{the max number of elements of a matrix we display
  as regular matrix with \code{image} or \code{display}. Larger matrices are represented
  as dots only.}

\item{trivalues}{a flag whether to return the structure (default,
  \code{FALSE} or the values themselves (\code{TRUE}) when returning the
  upper and lower triangular part of a matrix.}

\item{cex}{default dot size for  \code{image} or \code{display}.}

\item{version}{a list containing the (major and minor) version number
  and the release date.}

\item{bcksl=TRUE}{default parameter for "solve" routines. \code{FALSE}
  would solve the system without using the permutation.}
  
\item{safemode=TRUE}{verify double and integer formats when constructing
  \code{spam} objects. \code{TRUE} is safe but only slightly slower.}

}
\details{
  The user should only change \code{eps}, \code{printsize},
  \code{imagesize} and \code{drop}.\cr  
  \code{bcksl} and \code{safemode} are rather developing flags.
  It is the users responsibility to use a \code{safemode=FALSE}. In many
  cases, all runs smoothly but there might be instances (especially when
  assigning directly to slots) that crash \R. The verification is done
  with the functions \code{icheck} and \code{dcheck}.
}
%\references{
%}
\seealso{\code{\link{print}},
  \code{\link{display}},  \code{\link{image}}, \code{\link{upper.tri}}.
}
\examples{
smat <- diag.spam( 1:8)
smat
.Spam$printsize <- 49
smat


%# safemode is only slightly faster...
%tmp <- rep( rnorm(10^5), 100)
%.Spam$safemode <- FALSE
%system.time( dcheck( tmp))
%.Spam$safemode <- TRUE
%system.time( dcheck( tmp))


}
\author{Reinhard Furrer}
\keyword{IO}
