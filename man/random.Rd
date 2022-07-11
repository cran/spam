% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/random.Rd.                          %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\name{random}
\alias{spam_random}
\title{Create Random Sparse Matrices}
\description{Creates random spam matrix given the dimension and other parameters.}
\usage{
spam_random(nrow = 1L, ncol = nrow, density = 0.5, distribution = NULL, digits = NULL,
            sym = FALSE, spd = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{nrow}{integer value for the number of rows for the \code{spam} matrix to create.}
  \item{ncol}{integer value for the number of columns. The default value is the same as \code{nrow}.}
  \item{density}{A numeric value between 0 and 1 specifying the approximate density of matrix.
    If equal to zero the \code{spam} matrix contains only zeros and if equal to 1 the \code{spam} matrix is full.}
  \item{distribution}{a random number generating distribution function to sample the entries of the \code{spam} matrix.
    The function must have an argument with the name \code{n}, possible candidates are \code{rnorm}, \code{rexp}, \code{rpois}, \code{rweibull}, etc..}
  \item{...}{possible additional arguments for the distribution function if specified with \code{distribution}.}
  \item{digits}{an integer value for the number of digits the entries should be rounded.}
  \item{sym}{logical value to specify symmetry of the \code{spam} matrix.}
  \item{spd}{logical value to specify positive definitness of the \code{spam} matrix, via diagonal dominace criteria.
    Note, if \code{spd} TRUE, then \code{sym} is overwritten to \code{TRUE} in any case.}
  \item{verbose}{logical value to specify verbose statments of the function.}
}
\value{A random matrix in \code{spam} format.}
\details{To create a random spam64 matrix, set \code{options(spam.force64 = TRUE)}.}
%\references{}
\seealso{\code{\link{spam-class}} and \code{\link{display.spam}}}
\examples{
set.seed(42)
rspam <- spam_random(500, digits = 2, distribution = rnorm, sd = 2, mean = 10, density = .01)
display.spam(rspam, cex = 2)
}
\author{Florian Gerber, Roman Flury, Reinhard Furrer}

