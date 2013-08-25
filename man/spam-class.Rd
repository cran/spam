% This is file ../spam/man/spam-class.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{spam-class}
\docType{class}
\alias{spam.class}
\alias{spam-class}
\alias{as.matrix,spam-method}
\alias{as.matrix.spam}
\alias{[<-,spam,ANY,ANY,ANY-method}
\alias{[<-,spam,matrix,matrix,ANY-method}
\alias{[<-,spam,matrix,missing,ANY-method}
\alias{[<-,spam,missing,missing,ANY-method}
\alias{[<-,spam,missing,vector,ANY-method}
\alias{[<-,spam,missing,vector,spam-method}
\alias{[<-,spam,spam,missing,ANY-method}
\alias{[<-,spam,vector,missing,ANY-method}
\alias{[<-,spam,vector,missing,spam-method}
\alias{[<-,spam,vector,vector,ANY-method}
\alias{[<-,spam,vector,vector,spam-method}
\alias{[<-,spam,missing,missing,numeric-method}
\alias{[<-,spam,missing,vector,numeric-method}
\alias{[<-,spam,vector,missing,numeric-method}
\alias{[<-,spam,vector,vector,numeric-method}
\alias{[<-,spam,matrix,missing,numeric-method}
\alias{[<-,spam,matrix,matrix,numeric-method}
\alias{[<-,spam,spam,missing,numeric-method}
\alias{[<-,spam,ANY-method}
\alias{[,spam,missing,missing,ANY-method}
\alias{[,spam,missing,vector,ANY-method}
\alias{[,spam,vector,missing,logical-method}
\alias{[,spam,vector,missing,missing-method}
\alias{[,spam,vector,vector,ANY-method}
\alias{[,spam,matrix,missing,missing-method}
\alias{[,spam,matrix,missing,logical-method}
\alias{[,spam,matrix,matrix,ANY-method}
\alias{[,spam,spam,missing,ANY-method}
\alias{[,spam,ANY,ANY,ANY-method}
\alias{Arith,spam,missing-method}
\alias{Arith,spam,numeric-method}
\alias{Arith,numeric,spam-method}
\alias{Compare,numeric,spam-method}
\alias{Compare,spam,numeric-method}
\alias{c,spam-method}
\alias{t,spam-method}
\alias{length<-,spam-method}
\alias{length,spam-method}
\alias{lower.tri,spam-method}
\alias{kronecker,spam,ANY-method}
\alias{kronecker,ANY,spam-method}
\alias{kronecker,spam,spam-method}
\alias{plot,spam,missing-method}
\alias{plot,spam,spam-method}
\alias{show,spam-method}
\alias{upper.tri,spam-method}

\title{Class "spam"}
\description{The \code{spam} class is a representation of sparse matrices.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("spam", entries, colindices,
    rowpointes, dimension)}.
  The standard "old Yale sparse format" is used to store sparse matrices.\cr
  The matrix
\code{x} is stored in row form.  The first element of row \code{i} is
\code{x@rowpointers[i]}. The length of row \code{i} is determined by
\code{x@rowpointers[i+1]-x@rowpointers[i]}. The column indices of \code{x} are stored in
the \code{x@colindices} vector. The column index for element \code{x@entries[k]} is
\code{x@colindices[k]}.}
\section{Slots}{
  \describe{
    \item{\code{entries}:}{Object of class \code{"numeric"} contains the
      nonzero values. }
    \item{\code{colindices}:}{Object of class \code{"integer"} ordered indices
      of the nonzero values. }
    \item{\code{rowpointers}:}{Object of class \code{"integer"} pointer to the beginning
    of each row in the arrays \code{entries} and \code{colindices}.}
    \item{\code{dimension}:}{Object of class \code{"integer"} specifying
    the dimension of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{as.matrix}{\code{signature(x = "spam")}:
      transforming a sparse
      matrix into a regular matrix.}
    \item{as.spam}{\code{signature(x = "spam")}:
      cleaning of a sparse matrix.}
    \item{[<-}{\code{signature(x = "spam", i,j, value)}:
      assigning a
      sparse matrix. The negative vectors are not implemented yet. }
    \item{[}{\code{signature(x = "spam", i, j)}:
      subsetting a
      sparse matrix. The negative vectors are not implemented yet. }
    \item{\%*\%}{\code{signature(x, y)}:
      matrix multiplication, all combinations of sparse with full
      matrices or vectors are implemented.}
    \item{c}{\code{signature(x = "spam")}:
      vectorizes the sparse matrix and takes account of the zeros. Hence
    the lenght of the result is \code{prod(dim(x))}.}
    \item{cbind}{\code{signature(x = "spam")}: binds sparse matrices, see \code{\link{cbind}} for details.}
    \item{chol}{\code{signature(x = "spam")}:
    see \code{\link{chol}} for details.}
    \item{diag}{\code{signature(x = "spam")}: 
    see \code{\link{diag}} for details.}
  \item{dim<-}{\code{signature(x = "spam")}: truncates or augments
    the matrix see \code{\link{dim}} for details.}
    \item{dim}{\code{signature(x = "spam")}: gives the dimension of the
      sparse matrix.}
    \item{image}{\code{signature(x = "spam")}:
      see \code{\link{image}} for details.}
    \item{display}{\code{signature(x = "spam")}:
      see \code{\link{display}} for details.}
%    \item{initialize}{\code{signature(.Object = "spam")}: ... }
    \item{length<-}{\code{signature(x = "spam")}: Is not implemented and
      causes an error.}
    \item{length}{\code{signature(x = "spam")}: gives the number of
      non-zero elements.}
    \item{lower.tri}{\code{signature(x = "spam")}: see \code{\link{lower.tri}} for details. }
    \item{Math}{\code{signature(x = "spam")}: see \code{\link{Math}} for details.}
    \item{Math2}{\code{signature(x = "spam")}: see \code{\link{Math2}} for details. }
    \item{norm}{\code{signature(x = "spam")}: calculates the norm of a matrix.}
    \item{plot}{\code{signature(x = "spam", y)}: same functionality as
      the ordinary \code{plot}. }
    \item{print}{\code{signature(x = "spam")}: see \code{\link{print}} for details.}
    \item{rbind}{\code{signature(x = "spam")}: binds sparse
      matrices, see \code{\link{cbind}} for details.}
%    \item{show}{\code{signature(object = "spam")}: ... }
    \item{solve}{\code{signature(a = "spam")}: see \code{\link{solve}} for details.}
    \item{summary}{\code{signature(object = "spam")}: small summary
      statement of the sparse matrix.}
    \item{Summary}{\code{signature(x = "spam")}:
      All functions of the \code{Summary} class (like \code{min}, \code{max},
      \code{range}...) operate on the vector \code{x@entries} and return the
      result thereof. See Examples or \code{\link{Summary}} for details.
    }
    \item{t}{\code{signature(x = "spam")}: transpose of a sparse matrix.}
    \item{upper.tri}{\code{signature(x = "spam")}: see \code{\link{lower.tri}} for details.}
  }
}
\section{Details}{
  The compressed sparse row (CSR) format is often described with the
  vectors \code{a}, \code{ia}, \code{ja}. To be a bit more
  comprehensive, we have chosen longer slot names.
}
\section{Note}{The slots \code{colindices} and \code{rowpointers} are
  tested for proper integer assignments. This is not true for
  \code{entries}. 
}
%\references{}
\author{Reinhard Furrer, some of the Fortran code is based on A. George,
  J. Liu, E. S. Ng, B.W Peyton and Y. Saad (alphabetical)}
%\note{}
%\seealso{}
\examples{
showMethods("as.spam")


smat <- diag.spam(runif(15))
range(smat)
cos(smat)
}
\keyword{classes}
