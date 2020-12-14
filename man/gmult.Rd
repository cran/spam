% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/gmult.Rd.                           %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{gmult}
\alias{gmult}
\title{Generalized Multiplication}
\usage{
gmult(x, splits, fact)
}
\description{
Multiplying specific submatrices of a spam matrix with different factors.
}
\arguments{
  \item{x}{a spam matrix.}
  \item{splits}{vector of how to split the matrix into submatrices.
                It starts with \code{1} and ends with \code{max(dim(X))+1}.}
  \item{fact}{matrix of factors to multiply submatrices defined by splits.
              Dimensions of \code{fact} must correspond to the \code{length(splits)-1}.}
}
\value{Spam matrix, where each specified submatrix is multiplied with a factor.}

\examples{
x <- spam(1, 15, 15)
print(x, minimal = FALSE)

splits <- c(1,2,8,ncol(x)+1) # divide matrix into 9 submatrices
fact <- array(1:9, c(3,3))   # multiply each submatrix with a different factor
oF <- gmult(x, splits, fact)
print(oF, minimal = FALSE)
}

\author{Florian Gerber, Roman Flury}
\keyword{array}
\keyword{algebra}
