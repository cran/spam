\name{spam-class}
\docType{class}
\alias{spam.class}
\alias{spam-class}
\alias{as.matrix,spam-method}
\alias{as.matrix.spam}
\alias{[<-,spam,missing,missing,numeric-method}
\alias{[<-,spam,missing,vector,numeric-method}
\alias{[<-,spam,vector,missing,numeric-method}
\alias{[<-,spam,vector,vector,numeric-method}
\alias{[<-,spam,matrix,missing,numeric-method}
\alias{[<-,spam,matrix,matrix,numeric-method}
\alias{[<-,spam,spam,missing,numeric-method}
\alias{[<-,spam,ANY-method}
\alias{[,spam,missing,missing-method}
\alias{[,spam,missing,vector-method}
\alias{[,spam,vector,missing-method}
\alias{[,spam,vector,vector-method}
\alias{[,spam,matrix,missing-method}
\alias{[,spam,matrix,matrix-method}
\alias{[,spam,spam,missing-method}
\alias{[,spam,ANY-method}
\alias{Arith,spam,missing-method}
\alias{Arith,spam,numeric-method}
\alias{Arith,numeric,spam-method}
\alias{Compare,numeric,spam-method}
\alias{Compare,spam,numeric-method}
\alias{c,spam-method}
\alias{cbind,spam-method}
\alias{dim<-,spam-method}
\alias{image,spam-method}
\alias{length<-,spam-method}
\alias{length,spam-method}
\alias{lower.tri,spam-method}
\alias{plot,spam,missing-method}
\alias{plot,spam,spam-method}
\alias{rbind,spam-method}
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
      nonzero values }
    \item{\code{colindices}:}{Object of class \code{"integer"} ordered indices
      of the nonzero values }
    \item{\code{rowpointers}:}{Object of class \code{"integer"} pointer to the beginning
    of each row in the arrays \code{entries} and \code{colindices}}
    \item{\code{dimension}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{as.matrix}{\code{signature(x = "spam")}: ... }
    \item{as.spam}{\code{signature(x = "spam")}: ... }
    \item{[<-}{\code{signature(x = "spam", i = "missing", j = "missing", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "spam", i = "missing", j = "vector", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "spam", i = "vector", j = "missing", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "spam", i = "vector", j = "vector", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "spam", i = "matrix", j = "missing", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "spam", i = "matrix", j = "matrix", value = "numeric")}: ... }
    \item{[<-}{\code{signature(x = "spam", i = "spam", j = "missing", value = "numeric")}: ... }
    \item{[}{\code{signature(x = "spam", i = "missing", j = "missing")}: ... }
    \item{[}{\code{signature(x = "spam", i = "missing", j = "vector")}: ... }
    \item{[}{\code{signature(x = "spam", i = "vector", j = "missing")}: ... }
    \item{[}{\code{signature(x = "spam", i = "vector", j = "vector")}: ... }
    \item{[}{\code{signature(x = "spam", i = "matrix", j = "missing")}: ... }
    \item{[}{\code{signature(x = "spam", i = "matrix", j = "matrix")}: ... }
    \item{[}{\code{signature(x = "spam", i = "spam", j = "missing")}: ... }
    \item{[}{\code{signature(x = "spam", i = "ANY")}: ... }
    \item{\%*\%}{\code{signature(x = "spam", y = "spam")}: ... }
    \item{\%*\%}{\code{signature(x = "spam", y = "matrix")}: ... }
    \item{\%*\%}{\code{signature(x = "spam", y = "numeric")}: ... }
    \item{\%*\%}{\code{signature(x = "matrix", y = "spam")}: ... }
    \item{\%*\%}{\code{signature(x = "numeric", y = "spam")}: ... }
    \item{c}{\code{signature(x = "spam")}: ... }
    \item{cbind}{\code{signature(deparse.level = "spam")}: ... }
    \item{chol}{\code{signature(x = "spam")}: ... }
    \item{diag}{\code{signature(x = "spam")}: ... }
    \item{dim<-}{\code{signature(x = "spam")}: ... }
    \item{dim}{\code{signature(x = "spam")}: ... }
    \item{image}{\code{signature(x = "spam")}: ... }
    \item{initialize}{\code{signature(.Object = "spam")}: ... }
    \item{length<-}{\code{signature(x = "spam")}: ... }
    \item{length}{\code{signature(x = "spam")}: ... }
    \item{lower.tri}{\code{signature(x = "spam")}: ... }
    \item{Math}{\code{signature(x = "spam")}: ... }
    \item{norm}{\code{signature(x = "spam")}: ... }
    \item{plot}{\code{signature(x = "spam", y = "missing")}: ... }
    \item{plot}{\code{signature(x = "spam", y = "spam")}: ... }
    \item{print}{\code{signature(x = "spam")}: ... }
    \item{rbind}{\code{signature(deparse.level = "spam")}: ... }
    \item{show}{\code{signature(object = "spam")}: ... }
    \item{solve}{\code{signature(a = "spam")}: ... }
    \item{summary}{\code{signature(object = "spam")}: ... }
    \item{t}{\code{signature(x = "spam")}: ... }
    \item{upper.tri}{\code{signature(x = "spam")}: ... }
  }
}
%\references{}
\author{Reinhard Furrer}
%\note{}
%\seealso{}
%\examples{}
\keyword{classes}
