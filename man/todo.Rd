% This is file ../spam0.20-3/man/todo.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{todo}
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
  
  - extend demo(s)\cr
  
  - write vignette\cr

  - complete help files\cr
  %Especially, there are still a few functions
  %commented in various 'usage' sections (to avoid warnings when checking
  %the library).\cr

  - extend basic matrix operatation, comparisons, etc: \cr
    \code{unique}, \code{duplicated}, ... 
  
  - improve subsetting via row extraction, incorporate matrix permutation\cr % xtrows
  
%  - extend methods for \code{spam.chol.NgPeyton}: \cr
%  e.g., \code{diag}  => done
  
  - implement other Cholesky routines (one eye glances to the LDL library). \cr

  - what about an LU/SVD decomposition?\cr
  
%   - implement \code{backsolve} and \code{forwardsolve} for \code{spam} class.\cr  => done

%  - it might be useful to implement certain lines of code in Fortran:
%  say \code{rbind}, \code{cbind}, ...\cr  => done

  
  Any other items are welcome (\email{reinhard.furrer@math.uzh.ch}).\cr
}
\keyword{internal}
