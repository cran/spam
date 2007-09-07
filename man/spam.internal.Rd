\name{spam internal}
\alias{auxiliarycholNgPeyton}
\alias{dcheck}
\alias{icheck}
\title{
  Spam internal and auxiliary functions
}
\description{
  The functions listed below are auxiliary functions and are not
  exported by the NAMESPACE. The user should not
  require to call these directly. 
}
%\usage{
%auxiliarycholNgPeyton(x, nrow, memory)
%setdiagspam(x, value)
%}
\details{The functions are listed here for a better understanding of the
  code (to fulfill the tutorial style paradigm).\cr

  \code{auxiliarycholNgPeyton(x, nrow, memory)} The essential
  Cholesky factorization used by \code{chol}, \code{det}.

  \code{dcheck(x)}, \code{icheck(x)} testing and forcing of \code{x}
  to doubles and
  integers. These functions are used when calling Fortran routines.
  }
\keyword{internal}
