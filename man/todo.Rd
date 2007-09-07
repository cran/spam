\name{spam todo}
\alias{spam.todo}
\alias{todo.spam}
\alias{todo}
\title{Small "ToDo" list}
\description{List of what needs to be done within \code{spam}
}
\usage{
  todo()
}
\details{This is a non exhaustive list of where we need to work on spam
  (of course the list is in random order):\cr

  - nicer \code{.Spam} handling (without causing the warning, results
  from my poor NAMESPACE understanding)\cr
  
  - extend demo(s)\cr
  
  - write vignette\cr

  - complete help files. Especially, there are still a few functions
  commented in various 'usage' sections (to avoid warnings when checking
  the library).\cr

  - extend basic matrix operatation and comparisons, implement kronecker product \cr
  
  - implement preordering and symbolic factorization in \code{chol}\cr
  
  - implement other Cholesky routines (one eye glances to the LDL
  library). \cr
  % fill-in ratios...
  
  - implement \code{backsolve} and \code{forwardsolve} for \code{spam} class\cr

  Any other items are welcome (\email{rfurrer@mines.edu}).\cr

  
  On a lower level, it might be useful to implement certain lines of
  code in Fortran:\cr

  - Loop in \code{cbind.spam}\cr
  
  
}
\keyword{internal}
