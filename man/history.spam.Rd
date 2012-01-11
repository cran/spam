% This is file ../spam0.28-0/man/history.spam.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   


\name{history}
\alias{history.spam}
\alias{spam.history}
\alias{history}
\title{Brief Overview of the History}
\description{Brief overview of \code{spam}s  history}
\usage{
  spam.history()
}
\details{Many years ago, I required sparse matrix algebra capacities in
  \R for my research. The packages available at that time where not
  necessarily satisfying for my purpose. Hence, I started writing my own
  functions, somewhat inspired by the package \code{SparseM}.

  Over the years, I decided to write an \R package, called \code{spam}.
  In mean time, we have tree significant sparse matrix packages posted
  on CRAN. I believe that all three have their own niche of users.

  I believe that \code{spam} is stable, useful but definitely not fully
  matured. I appreciate all comments and remarks.


  There is no need to maintain here a complete list of changes. Older
  versions are available on CRAN (>0.11), each coming with a proper
  history file.
	
}
\references{See also a more detailed description at
  \url{http://www.math.uzh.ch/furrer/software/spam/history.shtml}.}
\seealso{\code{\link{todo}}.}
\examples{
\dontrun{
## A complete history of changes
file.show(system.file("NEWS", package = "spam"))
file.show(system.file("0NEWS", package = "spam"))
}
}
\keyword{internal}
