% This is file ../spam/man/options.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
   


\name{options}
\alias{options.spam}
\alias{spam.options}
\alias{spam.getOption}
\alias{.Spam}
\title{Options Settings}
\description{
  Allow the user to set and examine a variety of \emph{options}
  which affect the way in which \R computes and displays sparse
  matrix results.
}
\usage{
spam.options(\dots)

spam.getOption(x)

}
\arguments{
  \item{\dots}{any options can be defined, using \code{name = value} or
    by passing a list of such tagged values.  However, only the ones
    below are used in \code{spam}.\cr
    Further, \code{spam.options('name') == spam.options()['name']}, see the example.
  }
  \item{x}{a character string holding an option name.}
}
\details{
  Invoking \code{spam.options()} with no arguments returns a list with the
  current values of the options. To access the value of a single option, one should
  use \code{spam.getOption("eps")}, e.g., rather than
  \code{spam.options("eps")} which is a \emph{list} of length one.\cr

  Internally, the options are kept in the list \code{.Spam}.\cr

  Of course, printing is still subordinate to
  \code{getOption("max.print")} or similar options.
}
\value{
  For \code{spam.getOption}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.

  For \code{spam.options()}, a list of all set options sorted by category.  For
  \code{spam.options(name)}, a list of length one containing the set value,
  or \code{NULL} if it is unset.  For uses setting one or more options,
  a list with the previous values of the options changed (returned
  invisibly).
}
\section{Options used}{
  A short description with the default values follows.
  \describe{
\item{\code{eps=.Machine$double.eps}:}{values smaller than this are
  considered as zero. This is only used when creating spam objects.}

\item{\code{drop=FALSE}:}{default parameter for \code{drop} when subsetting}

\item{\code{printsize=100}:}{the max number of elements of a matrix which we
  display as regular matrix.}

\item{\code{imagesize=10000}:}{the max number of elements of a matrix we display
  as regular matrix with \code{image} or \code{display}. Larger matrices are represented
  as dots only.}

\item{\code{cex=1200}:}{default dot size for  \code{image} or \code{display}.}

\item{\code{structurebased=FALSE}:}{should operations be carried out on
  the nonzero entries (the structure) or including the zeros.}

\item{\code{inefficiencywarning=1e6}:}{issue a warning when inefficient
  operations are performed and the matrix exceeds the specified size.
  Valid value is a postive integer or a logical. \code{TRUE} corresponds
  to 1 (always), \code{FALSE} to \code{Inf}.}

\item{\code{trivalues=FALSE}:}{a flag whether to return the structure 
  (\code{FALSE}) or the values themselves (\code{TRUE}) when returning the
  upper and lower triangular part of a matrix.}

\item{\code{listmethod="PE"}:}{algorithm for \code{spam.list}. Default
  is suggestion by Paul Eilers (thanks). Any other specification uses a
  bubble sort algorithm which is only slightly faster for very sparse matrices. }

\item{\code{dopivoting=TRUE}:}{default parameter for "solve" routines. \code{FALSE}
  would solve the system without using the permutation.}
  
\item{\code{NAOK=FALSE}:}{logical determines if NA, NaN and Inf are
  allowed to Fortan. Setting to TRUE allows to work with these but
  full functionality has not been tested.}

\item{\code{safemodevalidity=TRUE}:}{logical determines if sanity check
  is peformed when constructing sparse matrices.  
  Default is safer but somewhat slower.}

\item{\code{cholsymmetrycheck=TRUE}:}{for the Cholesky factorization,
  verify if the matrix is symmetric.}
\item{\code{cholpivotcheck=TRUE}:}{for the Cholesky factorization,
  when passing a permutation, should a minimum set of checks be performed?}
\item{\code{cholupdatesingular="warning"}:}{for a Cholesky update, what
  happens if the matrix is singular: \code{"warning"} only and
  returning the not updated factor, \code{"error"}  or  return simply \code{"NULL"}.}
\item{\code{cholincreasefactor=c(1.25,1.25)}:}{If not enought memory
  could be allocated, these are the steps to increase it.}

\item{\code{nnznearestdistnnz=c(400^2,400)}:}{Memory allocation
  parameters for \code{nearest.dist}.}
\item{\code{nearestdistincreasefactor=1.25}:}{If not enought memory
  could be allocated, this is the step to increase it. }
}
}
\seealso{Functions influenced by these options include: \code{\link{print.spam}},
  \code{\link{display.spam}},  \code{\link{image.spam}}, \code{\link{upper.tri.spam}},
  \code{\link{chol.spam}}, \code{\link{nearest.dist}}, etc.\cr
  \code{\link{powerboost}}\cr
}
\author{\code{spam.options} is essentially identical to \code{sm.options}.}
\examples{
op <- spam.options()

# two ways of representing the options nicely.
utils::str(op)
noquote(format(op) )


smat <- diag.spam( 1:8)
smat
spam.options( printsize=49)
smat


# Reset to default values:
spam.options( eps=.Machine$double.eps, drop=FALSE,             
   printsize=100, imagesize=10000, cex=1200, 
   structurebased=FALSE, inefficiencywarning=1e6,
   trivalues=FALSE, listmethod="PE", 
   NAOK=FALSE, safemodevalidity=TRUE, dopivoting=TRUE, 
   cholsymmetrycheck=TRUE, cholpivotcheck=TRUE,
   cholupdatesingular="warning", cholincreasefactor=c(1.25,1.25),
   nearestdistincreasefactor=1.25, nearestdistnnz=c(400^2,400)     )

}
\keyword{IO}
\keyword{environment}
\keyword{error}
\keyword{print}
