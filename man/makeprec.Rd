% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/makeprec.Rd.                        %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{makeprec}
\alias{precmat.GMRFreglat}
\title{Create Precision Matrices}
\description{Creates precision matrices for gridded GMRF.}
\usage{
precmat.GMRFreglat(n,m, par, model = "m1p1",  eps = getOption("spam.eps"))
}
\arguments{
  \item{n}{first dimension of the grid.}
  \item{m}{second dimension of the grid.}
  \item{par}{parameters used to construct the matrix.}
  \item{model}{see details and examples.}
  \item{eps}{A tolerance parameter:  elements of \code{x} such that \code{abs(x)
    <= eps} set to zero. Defaults to \code{eps = getOption("spam.eps")}} 
}
\value{A \code{spam} matrix of dimension \code{prod(dims)}x\code{prod(dims)}.}
\details{The function should be illustrative on how to create precision
  matrices for gridded GMRF. Hence, no testing (positive definiteness is
  done).

  The model specification \code{"m"} determines the complexity and
  \code{"p"} the number of parameters. 
  
  Please see the examples on the meaning of the different models. 
    }
%\references{}
\seealso{\code{\link{precmat}}, \code{\link{toeplitz.spam}}, \code{\link{kronecker.spam}}}
\examples{
as.matrix(precmat.GMRFreglat(4, 3, c(.4),         'm1p1'))
as.matrix(precmat.GMRFreglat(4, 3, c(.4,.3),      'm1p2'))
as.matrix(precmat.GMRFreglat(4, 3, c(.4,.3,.2),   'm2p3'))
as.matrix(precmat.GMRFreglat(4, 3, c(.4,.3,.2,.1),'m2p4'))

# up to the diagonal, the following are equivalent:
cleanup( precmat.IGMRFreglat(3,4) -
             precmat.GMRFreglat(3,4,1, 'm1p1'))

}
\author{Reinhard Furrer}
\keyword{array}
\keyword{algebra}
