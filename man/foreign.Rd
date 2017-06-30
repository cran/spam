% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/foreign.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{foreign}
\alias{foreign}
\alias{as.spam.matrix.csr}
%\alias{as.matrix.csr.spam}
\alias{as.dgRMatrix.spam}
\alias{as.dgCMatrix.spam}
\alias{as.spam.dgRMatrix}
\alias{as.spam.dgCMatrix}
\title{Transformation to other sparse formats}
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
\seealso{ \code{\link{triplet}}, \code{\link[Matrix]{Matrix}} or   \code{\link[SparseM]{matrix.csr}}}
\examples{

\dontrun{
S <- diag.spam(4)
R <- as.dgRMatrix.spam( S) 
C <- as.dgCMatrix.spam( S)
as.spam.dgCMatrix(C) 
slotNames(U)
slotNames(R)
# For column oriented sparse formats a transpose does not the job,
# as the slot names change.


# as.spam(R) does not work.
}

\dontrun{
# for transformations between SparseM and spam:
as.matrix.csr.spam <- function(x,...) {
    newx <- new("matrix.csr")
    slot(newx,"ra",check=FALSE) <- x@entries
    slot(newx,"ja",check=FALSE) <- x@colindices
    slot(newx,"ia",check=FALSE) <- x@rowpointers
    slot(newx,"dimension",check=FALSE) <- x@dimension
    return(newx)
  }       
U <- as.matrix.csr.spam( S) 
}

\dontrun{
# a dataset contained in Matrix
data(KNex)
as.spam.dgCMatrix(KNex$mm) 

}
}
\author{Reinhard Furrer}
\keyword{array}
\keyword{manip}
