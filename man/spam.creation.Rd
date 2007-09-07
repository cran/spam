\name{spam}
\alias{spam.creation}
\alias{initialize,spam-method}
\alias{spam}
\alias{spam.numeric}
\alias{spam.spam}
\alias{spam,spam-method}
\alias{spam,numeric-method}
\alias{as.spam,matrix-method}
\alias{as.spam,numeric-method}
\alias{as.spam,spam-method}
\alias{as.spam,spam.chol.NgPeyton-method}
\alias{as.spam}
\alias{as.spam.spam}
\alias{as.spam.numeric}
\alias{as.spam.matrix}
\alias{as.spam.chol.NgPeyton}
\alias{as.spam.}
\alias{is.spam}
\alias{validspamobject}
\title{Sparse Matrix Class}
\description{
This group of functions evaluates and coerces changes in class structure.
}
\usage{
spam(x, nrow = 1, ncol = 1, eps = .Spam$eps)

as.spam(x, eps = .Spam$eps)

is.spam(x)
}
\arguments{
\item{x}{is a matrix, or vector object,  of either dense or sparse form}
\item{nrow}{number of rows of matrix }
\item{ncol}{number of columns of matrix }
\item{eps}{A tolerance parameter:  elements of \code{x} such that \code{abs(x) < eps} set to zero.
This argument is only relevant when coercing matrices from dense to sparse form. Defaults to
\code{eps = .Spam$eps} }

}
\details{
The function \code{as.spam} acts like \code{matrix} to coerce a vector object to
a sparse matrix object of class \code{spam}.  
}
%\note{}
%\references{}
\seealso{
\code{\link{SPAM}} general overview of the package.
}
\examples{
# old message, do not loop, when you create a large sparse matrix
nz <- 128
ln <- nz^2
smat <- spam(0,ln,ln)
system.time(for (i in 1:nz) smat[sample(ln,1),sample(ln,1)] <- i)
system.time(smat[cbind(sample(ln,nz),sample(ln,nz))] <- 1:nz)

getClass("spam")

}
\keyword{algebra}
