\name{SPAM}
\alias{overview}
\alias{SPAM}
\alias{Spam}
\title{SPArse Matrix Package}
\description{\code{spam} is a collection of functions for sparse matrix
  algebra. 
}
\section{Gereral overview}{What is spam and what is it not\cr

While \code{Matrix} seems an overshoot of classes and \code{SparseM}
focuses mainly on regression type problem, we provide a minimal set of
sparse matrix functions fully functional for everyday spatial statistics
life.\cr

Emphasis is given on a comprehensiv, simple, tutorial structure of the
code. The code is S4 based but (in a tutorial spirit) the functions are
in a S3 structure visible to the user (exported via NAMESPACE).\cr

Sparseness is used when handling large matrices. Hence, care has been
used to provide efficient and fast routines. Essentially, the functions
do not transform the sparse structure into full matrices to use standard
(available) functionality, followed by a back transform. We agree, more
operators, functions, etc. should eventually be implemented.


The spatial package \code{fields} uses \code{spam} as a required
package.
\cr
  
}
\references{\url{www.mines.edu/~rfurrer/software/spam/}
}
\seealso{See \code{\link{spam.class}} for a detailed class description,
  \code{\link{spam}} and \code{\link{spam.ops}} for creation,
  coercion and algebraic operations.\cr
  \code{demo(package='spam')} lists available demos.\cr 
  Related packages are \code{\link[fields]{fields}},
  \code{\link[Matrix]{Matrix}} and
  \code{\link[SparseM]{SparseM.ontology}}.
}
\examples{
}
\author{Reinhard Furrer}
\keyword{documentation}