% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/det.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

\name{det}
\alias{det,spam-method}
\alias{det,spam.chol.NgPeyton-method}
\alias{det.spam}
\alias{determinant}
%\alias{determinant,spam.chol.NgPeyton-method}
\alias{determinant,spam-method}
\alias{determinant.spam}
\alias{determinant.spam.chol}
\alias{determinant.spam.chol.NgPeyton}
\title{Calculate the determinant of a positive definite Sparse Matrix}
\description{\code{det} and \code{determinant} calculate the determinant of a
  positive definite sparse matrix. \code{determinant} returns separately the
  modulus of the determinant, optionally on the logarithm scale, and
  the sign of the determinant.

}
\usage{
# det(x, ...)
determinant(x, logarithm = TRUE, ...)
}
\arguments{
  \item{x}{sparse matrix of class \code{spam} or a Cholesky factor  of
    class \code{spam.chol.NgPeyton}.}
  \item{logarithm}{logical; if \code{TRUE} (default) return the logarithm of the
    modulus of the determinant.}
  \item{...}{Optional arguments. Examples include \code{method} argument
    and additional parameters used by the method.}
}
\value{For \code{det}, the determinant of \code{x}.  For \code{determinant}, a list with
     components

 \item{modulus}{a numeric value.  The modulus (absolute value) of the
          determinant if \code{logarithm} is \code{FALSE}; otherwise the
          logarithm of the modulus.}

   \item{sign}{integer; either +1 or -1 according to whether the determinant
          is positive or negative.}

}
\details{If the matrix is not positive definite, the function issues a
  warning and returns \code{NA}.

  The determinant is based on the product of the diagonal entries of a
  Cholesky factor, i.e. internally, a Cholesky decomposition is
  performed. By default, the NgPeyton algorithm with minimal degree
  ordering us used. To change the methods or supply additonal parameters
  to the Cholesky factorization function, see the help for \code{\link{chol}}.

  The determinant of a Cholesky factor is also defined.
}
\references{
Ng, E. G. and B. W. Peyton (1993) Block sparse Cholesky algorithms
on advanced uniprocessor computers, \emph{SIAM J. Sci. Comput.}, \bold{14},
1034--1056.
}
\seealso{\code{\link{chol.spam}}
}
\examples{
x <- spam( c(4,3,0,3,5,1,0,1,4),3)
det( x)
determinant( x)

det( chol( x))
}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
