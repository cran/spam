% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/xybind.Rd.                          %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

\name{cbind}
\alias{cbind.spam}
\alias{rbind.spam}
\alias{cbind}
\alias{rbind}
\alias{cbind,spam-method}
\alias{rbind,spam-method}
\title{Combine Sparse Matrices by Rows or Columns}
\description{Take a sequence of vector, matrix or \code{spam} object arguments and
     combine by \emph{c}olumns or \emph{r}ows, respectively.}
\usage{
# cbind(\dots, force64 = getOption("spam.force64"), deparse.level = 0)
# rbind(\dots, deparse.level = 0)
}
\arguments{
  \item{...}{vectors, matrices or \code{spam} objects. See \sQuote{Details} and \sQuote{Value}}
  \item{force64}{logical vector of length 1. If \code{TRUE}, a 64-bit
    spam matrix is returned in any case. If \code{FALSE}, a 32-bit
    matrix is returned when possible. }
  \item{deparse.level}{for compatibility reason here. Only \code{0} is implemented.}
}
\value{a \code{spam} object combining the \code{\dots} arguments
  column-wise or row-wise.  (Exception: if there are no inputs or all
  the inputs are \code{NULL}, the value is \code{NULL}.)}
\details{\code{rbind} and \code{cbind} are not exactly symmetric in how
  the objects are processed. The former is essentially an concatenation
  of the slots due to the sparse storage format. Different types of
  inputs are handled differently. The latter calls a
  Fortran routine after the input has been coerced to \code{spam}
  objects. \cr

  Only two objects at a time are processed. If more than two are
  present, a loop concatenates them successively. \cr

  A method is defined for a  \code{spam} object as first argument. 
    }
%\references{}
%\seealso{\code{\link{cbind,spam-method}}.}
\examples{
x <- cbind.spam(1:5,6)

y <- cbind(x, 7)

rbind( x, x)
# for some large matrices   t( cbind( t(x), t(x)))
# might be slightly faster:


}
\author{Reinhard Furrer}
\keyword{array}
\keyword{manip}
