% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/foreign.Rd.                         %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\name{foreign}
\alias{foreign}
\alias{as.spam.matrix.csr}
%\alias{as.matrix.csr.spam}
\alias{as.dgRMatrix.spam}
\alias{as.dgCMatrix.spam}
\alias{as.spam.dgRMatrix}
\alias{as.spam.dgCMatrix}
\title{Transformation to Other Sparse Formats}
\description{Transform between the \code{spam} sparse format to the
  \code{matrix.csr} format of \code{SparseM} and  \code{dgRMatrix} format of \code{Matrix}}
\usage{
as.spam.matrix.csr(x)
as.dgRMatrix.spam(x)
as.dgCMatrix.spam(x)
as.spam.dgRMatrix(x)
as.spam.dgCMatrix(x)
}
\arguments{
  \item{x}{sparse matrix of class \code{spam}, \code{matrix.csr},
    \code{dgRMatrix} or \code{dgCMatrix}.}
}
\value{According to the call, a sparse matrix of class  \code{spam},
  \code{matrix.csr}, \code{dgRMatrix} or \code{dgCMatrix}.}
\details{
  We do not provide any \code{S4} methods and because of the existing
  mechanism a standard \code{S3} does not work.\cr

  The functions are based on \code{require}.\cr

  Notice that \code{as.matrix.csr.spam} should read as \code{as."matrix.csr".spam}.
}
%\references{}
\seealso{\code{\link{triplet}}, \code{\link[Matrix]{Matrix}} or \code{matrix.csr} from package \pkg{SparseM}.}
\examples{

\dontrun{
S <- diag.spam(4)
R <- as.dgRMatrix.spam( S)
C <- as.dgCMatrix.spam( S)
as.spam.dgCMatrix(C)
slotNames(C)
slotNames(R)
# For column oriented sparse formats a transpose does not the job,
# as the slot names change.


# as.spam(R) does not work.
}

\dontrun{
# for transformations between SparseM and spam:

as.matrix.csr.spam <- function(x,...) {
  if (new("matrix.csr")) {
    newx <- new("matrix.csr")
    slot(newx,"ra",check=FALSE) <- x@entries
    slot(newx,"ja",check=FALSE) <- x@colindices
    slot(newx,"ia",check=FALSE) <- x@rowpointers
    slot(newx,"dimension",check=FALSE) <- x@dimension
    return(newx)
    }
  stop("function requires 'SparseM' package")
}
# then with `SparseM`:  as.matrix.csr.spam( spamobject )
}

\dontrun{
# a dataset contained in Matrix
data(KNex, package='Matrix')
summary( KN <- as.spam.dgCMatrix(KNex$mm) )

}
}
\author{Reinhard Furrer}
\keyword{array}
\keyword{manip}
