% This is file spam0.15-3/man/options.Rd
% This file is part of the spam package, 
%      http://www.mines.edu/~rfurrer/software/spam/
% written and maintained by Reinhard Furrer.





\name{options}
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
    below are used in \code{spam}.
    Further, \code{spam.options('name') == spam.options()['name']}, see the example.
  }
  \item{x}{a character string holding an option name.}
}
\details{
  Invoking \code{spam.options()} with no arguments returns a list with the
  current values of the options. To access the value of a single option, one should
  use \code{spam.getOption("eps")}, e.g., rather than
  \code{spam.options("eps")} which is a \emph{list} of length one.\cr

  Internally, the options are kept in the list \code{.Spam}.
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

\item{\code{trivalues=FALSE}:}{a flag whether to return the structure 
  (\code{FALSE}) or the values themselves (\code{TRUE}) when returning the
  upper and lower triangular part of a matrix.}

\item{\code{cex=1200}:}{default dot size for  \code{image} or \code{display}.}

\item{\code{dopivoting=TRUE}:}{default parameter for "solve" routines. \code{FALSE}
  would solve the system without using the permutation.}
  
\item{\code{safemode=c(TRUE,TRUE,TRUE)}:}{The logicals are determine (1) 
  verify double and integer formats when constructing
  \code{spam} objects (2) quick sanity check when constructing sparse
  matrices (3) testing for NAs in Fortan calls.
  \code{TRUE}s are safer but slightly slower. The most relevant speedup
  is the second flag.}

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
\seealso{\code{\link{print}},
  \code{\link{display}},  \code{\link{image}}, \code{\link{upper.tri}},
  \code{\link{chol}}, \code{\link{nearest.dist}}, etc.
}
\author{\code{spam.options} is essentially identical to \code{sm.options}.}
\examples{
op <- spam.options()

# two ways of representing the options nicely.
utils::str(op)
noquote(unlist(format(op)) )


smat <- diag.spam( 1:8)
smat
spam.options( printsize=49)
smat


# Reset to default values:
spam.options( eps=.Machine$double.eps, drop=FALSE,             
   printsize=100,  imagesize=10000,   cex=1200, 
   trivalues=FALSE,   safemode=c(TRUE,TRUE,TRUE), 
   dopivoting=TRUE,       cholsymmetrycheck=TRUE,  
   cholpivotcheck=TRUE, cholupdatesingular="warning",
   cholincreasefactor=c(1.25,1.25),
   nearestdistincreasefactor=1.25,
   nearestdistnnz=c(400^2,400)     )

}
\keyword{IO}
\keyword{environment}
\keyword{error}
\keyword{print}
