% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/options.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\name{options}
\alias{spam.options}
\alias{options.spam}
\alias{getOption.spam}
\title{Options Settings}
\description{
  Allow the user to set and examine a variety of \emph{options}
  which affect the way in which \R computes and displays sparse
  matrix results.
}
\details{
  Invoking \code{options()} with no arguments returns a list with the
  current values of the options. To access the value of a single option, one should
  use \code{getOption("spam.eps")}, e.g., rather than
  \code{options("spam.eps")} which is a \emph{list} of length one.\cr

  Of course, printing is still subordinate to
  \code{getOption("max.print")} or similar options.
}
\value{
  For \code{getOption}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.

  For \code{options()}, a list of all set options sorted by category.  For
  \code{options(name)}, a list of length one containing the set value,
  or \code{NULL} if it is unset.  For uses setting one or more options,
  a list with the previous values of the options changed (returned
  invisibly).
}
\section{Options used for the package \code{spam}}{
  A short description with the default values follows.
  \describe{
\item{\code{spam.eps=.Machine$double.eps}:}{values smaller than this are
  considered as zero. This is only used when creating spam objects.}

\item{\code{spam.drop=FALSE}:}{default parameter for \code{drop} when subsetting}

\item{\code{spam.printsize=100}:}{the max number of elements of a matrix which we
  display as regular matrix.}

\item{\code{spam.imagesize=10000}:}{the max number of elements of a matrix we display
  as regular matrix with \code{image} or \code{display}. Larger matrices are represented
  as dots only.}

\item{\code{spam.cex=1200}:}{default dot size for  \code{image} or \code{display}.}

\item{\code{spam.structurebased=FALSE}:}{should operations be carried out on
  the nonzero entries (the structure) or including the zeros.}

\item{\code{spam.inefficiencywarning=1e6}:}{issue a warning when inefficient
  operations are performed and the matrix exceeds the specified size.
  Valid value is a postive integer or a logical. \code{TRUE} corresponds
  to 1 (always), \code{FALSE} to \code{Inf}.}

\item{\code{spam.trivalues=FALSE}:}{a flag whether to return the structure
  (\code{FALSE}) or the values themselves (\code{TRUE}) when returning the
  upper and lower triangular part of a matrix.}

\item{\code{spam.listmethod="PE"}:}{algorithm for \code{spam.list}. Default
  is suggestion by Paul Eilers (thanks). Any other specification uses a
  bubble sort algorithm which is only slightly faster for very sparse matrices. }

\item{\code{spam.dopivoting=TRUE}:}{default parameter for "\code{solve}" routines. \code{FALSE}
  would solve the system without using the permutation.}

\item{\code{spam.NAOK=FALSE}:}{logical determines if \code{NA}, \code{NaN} and \code{Inf} are
  allowed to Fortan. Setting to \code{TRUE} allows to work with these but
  full functionality has not been tested.}

\item{\code{spam.safemodevalidity=TRUE}:}{logical determines if sanity check
  is peformed when constructing sparse matrices.
  Default is safer but somewhat slower.}

\item{\code{spam.cholsymmetrycheck=TRUE}:}{for the Cholesky factorization,
  verify if the matrix is symmetric.}
\item{\code{spam.cholpivotcheck=TRUE}:}{for the Cholesky factorization,
  when passing a permutation, should a minimum set of checks be performed?}
\item{\code{spam.cholupdatesingular="warning"}:}{for a Cholesky update, what
  happens if the matrix is singular: \code{"warning"} only and
  returning the not updated factor, \code{"error"}  or  return simply \code{"NULL"}.}
\item{\code{spam.cholincreasefactor=c(1.25,1.25)}:}{If not enought memory
  could be allocated, these are the steps to increase it.}

\item{\code{spam.nnznearestdistnnz=c(400^2,400)}:}{Memory allocation
  parameters for \code{nearest.dist}.}
\item{\code{spam.nearestdistincreasefactor=1.25}:}{If not enought memory
  could be allocated, this is the step to increase it. }
}
}
\seealso{Functions influenced by these options include: \code{\link{print.spam}},
  \code{\link{display.spam}},  \code{\link{image.spam}}, \code{\link{upper.tri.spam}},
  \code{\link{chol.spam}}, \code{\link{nearest.dist}}, etc.\cr
  \code{\link{powerboost}}\cr
}

\examples{
smat <- diag.spam(1:8)
smat
options(spam.printsize=49)
smat

# All spam options:
op <- options()[grep("spam", names(options()))]

# Reset to default values:
options(spam.eps=.Machine$double.eps,
        spam.drop=FALSE,
        spam.printsize=100,
        spam.imagesize=10000,
        spam.cex=1200,
        spam.structurebased=FALSE,
        spam.inefficiencywarning=1e6,
        spam.trivalues=FALSE,
        spam.listmethod="PE",
        spam.NAOK=FALSE,
        spam.safemodevalidity=TRUE,
        spam.dopivoting=TRUE,
        spam.cholsymmetrycheck=TRUE,
        spam.cholpivotcheck=TRUE,
        spam.cholupdatesingular="warning",
        spam.cholincreasefactor=c(1.25,1.25),
        spam.nearestdistincreasefactor=1.25,
        spam.nearestdistnnz=c(400^2,400))

}
\keyword{IO}
\keyword{environment}
\keyword{error}
\keyword{print}
