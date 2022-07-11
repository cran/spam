% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/eigen.Rd.                           %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{eigen}
\alias{eigen.spam}
\alias{eigen_approx}
\title{Eigenvalues for Sparse Matrices}
\description{
  Functions to calculate eigenvalues and eigenvectors of \code{sparse} matrices.
  It uses the value of \code{spam.options("inefficiencywarning")} to dispatch between \code{base::eigen()} or the Implicitly Restarted Arnoldi Process, using 'ARPACK'.

  \code{eigen.spam} is a wrapper function of \code{eigen_approx} and transforms its output to \code{base::eigen} like.
}
\note{
The user is advised to choose the \code{control} options carefully, see \sQuote{Details} for more information.
}
\usage{
eigen.spam(x, nev = 10, symmetric, only.values = FALSE, control = list())
eigen_approx(x, nev, ncv, nitr, mode, only.values = FALSE, verbose = FALSE, f_routine)
}
\arguments{
\item{x}{a matrix of class \code{spam} whose \code{nev} eigenvalues and eigenvectors are to be computed.}
\item{nev}{number of eigenvalues to calculate.}
\item{symmetric}{if TRUE, the matrix is assumed to be symmetric.}
\item{only.values}{if TRUE, only \code{nev} eigenvalues are computed and returned, otherwise \code{nev} eigenvalues and eigenvectors are returned.}
\item{control}{additional options, see \sQuote{Details}.}
\item{ncv}{see \sQuote{Details}, use the \code{control} option for \code{eigen.spam}.}
\item{nitr}{see \sQuote{Details}, use the \code{control} option for \code{eigen.spam}.}
\item{mode}{see \sQuote{Details}, use the \code{control} option for \code{eigen.spam}.}
\item{verbose}{see \sQuote{Details}, use the \code{control} option for \code{eigen.spam}.}
\item{f_routine}{only for \code{eigen_approx}, to call the Fortran routine for symmetric matrices set this option to "ds_eigen_f" and for non symmetric to "dn_eigen_f".}
}
\value{
  A vector of the length corresponding to the dimension of the input matrix.
  Containing the required \code{nev} eigenvalues.
  If requested also the corresponding eigenvectors.
  In the non symmetric case, the eigenvalues are returned in a matrix with a column containing the real parts and a column containing the imaginary parts of the eigenvalues.
  The eigenvectors are then returned in two matrices.}

\details{
  \describe{
    \item{\code{mode = "  "}:}{
      there are different modes available for this function, each mode returns a different range of eigenvalues.
      Also the available modes are dependent, whether the input matrix is symmetric or not:
      \describe{
        \item{\code{"LM"}:}{Eigenvalues with largest magnitude (sym, non sym), that is, largest eigenvalues in the Euclidean norm of complex numbers.}
        \item{\code{"SM"}:}{Eigenvalues with smallest magnitude (sym, non sym), that is, smallest eigenvalues in the Euclidean norm of complex numbers.}
        \item{\code{"LR"}:}{Eigenvalues with largest real part (non sym).}
        \item{\code{"SR"}:}{Eigenvalues with smallest real part (non sym).}
        \item{\code{"LI"}:}{Eigenvalues with largest imaginary part (non sym).}
        \item{\code{"SI"}:}{Eigenvalues with smallest imaginary part (non sym).}
        \item{\code{"LA"}:}{Eigenvalues with largest algebraic value (sym), that is, largest eigenvalues inclusive of any negative sign.}
        \item{\code{"SA"}:}{Eigenvalues with smallest algebraic value (syn), that is, smallest eigenvalues inclusive of any negative sign.}
      }% describe
    }% item
    \item{\code{ncv}:}{
      the largest number of basis vectors that will be used in the Implicitly Restarted Arnoldi Process.
      Work per major iteration is proportional to x@dimension[1]*ncv*ncv.
      The default is set if \code{symmetric} to min(x@dimension[1] + 1, max(2 * nev + 1, 200)) or else to min(x@dimension[1] - 1, max(2 * nev + 1, 100)).
      Note, this value should not be chosen arbitrary large, but slightly larger than \code{nev}.
      Otherwise it could lead to memory allocation problems.}
    \item{\code{nitr}:}{
      the maximum number of iterations.
      The default is set to \code{ncv + 1000}}
    \item{\code{spamflag = FALSE}:}{
      if TRUE, the Implicitly Restarted Arnoldi Process is used,
  independent of the dimension of the respective matrix (provided matrix
  is larger than 10x10).}
    \item{\code{verbose = FALSE}:}{
      print additional information.}
    \item{\code{cmplxeps}:}{
      threshold to determine whether a double value is zero, while transforming the ARPACK output to R class complex.
      The default is set to \code{.Machine$double.eps}.}
}% describe
}% details


\references{Lehoucq, R. B. and Sorensen, D. C. and Yang, C. (1997) \emph{ARPACK Users Guide: Solution of Large Scale Eigenvalue Problems by Implicitly Restarted Arnoldi Methods}.}
\seealso{Option \code{"inefficiencywarning"} in \code{\link{spam.options}} and \code{\link{spam_random}}.
}
\examples{
set.seed(81)
rspam <- spam_random(42^2, density = .0001, spd = TRUE)
SPD <- eigen.spam(rspam, nev = 18, control = list(mode = "SM"),
                  only.values = TRUE)

any(SPD$values <= 0, na.rm = TRUE)
isSymmetric(rspam)
# hence the matrix is symmetric positiv definit

rspam2 <- spam_random(50^2, density = .0001, spd = FALSE, sym = TRUE,
                      distribution = rpois, lambda = 2)
SNPD <- eigen.spam(rspam2, nev = 18, control = list(mode = "SM"),
                    only.values = TRUE)

any(SNPD$values <= 0, na.rm = TRUE)
isSymmetric(rspam2)
# hence the matrix is symmetric but not positiv definit
}
\author{Roman Flury, Reinhard Furrer}
\keyword{algebra}

