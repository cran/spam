% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/bandwidth.Rd.                       %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   










\name{bandwidth}
\alias{bandwidth}
\title{Bandwidth of a Sparse Matrix}
\description{Returns the lower and upper bandwidth of a sparse matrix}
\usage{
bandwidth(A)
}
\arguments{
  \item{A}{spam object}
}
\details{The matrix does not need to be diagonal. Values can be negative
  indicating the the matrix contains a band cinfined in the
  upper or lower triangular part.
	}
\value{Integer vector containing the lower and upper bandwidth}
%\references{}
\seealso{
	\code{\link{diag.spam}}.
}
\examples{
bandwidth(spam(c(0, 1), 3, 2))

bandwidth(spam(c(0, 0, 1, rep(0, 9)), 4, 3))

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
