% This is file ../spam0.27-0/man/makeprec.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   



\name{makeprec}
\alias{make.prec}
\title{Create Precision Matrices}
\description{Creates precision matrices for gridded GMRF.}
\usage{
make.prec(par, dims, model = "m1p1",  eps = .Spam$eps)
}
\arguments{
  \item{par}{parameters used to construct the matrix.}
  \item{dims}{dimension of the grid.}
  \item{model}{see details and examples.}
  \item{eps}{A tolerance parameter:  elements of 'x' such that 'abs(x)
    <= eps' set to zero. Defaults to 'eps = .Spam$eps'} 
}
\value{A \code{spam} matrix of dimension \code{prod(dims)}x\code{prod(dims)}.}
\details{The function should be illustrative on how to create precision
  matrices for gridded GMRF. Hence, no testing (positive definiteness is
  done). Please see the examples on the meaning of the different models. 
    }
%\references{}
\seealso{\code{\link{toeplitz.spam}}, \code{\link{kronecker.spam}}}
\examples{
as.matrix(make.prec(c(.4),         c(4,3),'m1p1'))
as.matrix(make.prec(c(.4,.3),      c(4,3),'m1p2'))
as.matrix(make.prec(c(.4,.3,.2),   c(4,3),'m2p3'))
as.matrix(make.prec(c(.4,.3,.2,.1),c(4,3),'m2p4'))

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
