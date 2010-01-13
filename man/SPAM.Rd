% This is file ../spam0.20-3/man/SPAM.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   










\name{. SPAM .}
\alias{overview}
\alias{SPAM}
\alias{Spam}
\docType{package}
\title{SPArse Matrix Package}
\description{\code{spam} is a collection of functions for sparse matrix
  algebra. 
}
\section{Gereral overview}{What is spam and what is it not:\cr

While \code{Matrix} seems an overshoot of classes and \code{SparseM}
focuses mainly on regression type problem, we provide a minimal set of
sparse matrix functions fully functional for everyday spatial statistics
life. There is however some emphasize on Markov chain Monte Carlo type
calculations within the framework of (Gaussian) Markov random fields.
\cr

Emphasis is given on a comprehensive, simple, tutorial structure of the
code. The code is S4 based but (in a tutorial spirit) the functions are
in a S3 structure visible to the user (exported via \code{NAMESPACE}).\cr

There exist many methods for sparse matrices that work identically as in
the case of ordinary matrices. All the methods are discussed in the help
and can be accessed directly via a \code{*.spam} concatenation to the
function. For example, \code{help{chol.spam}} calls the help directly,
whereas with  \code{help{chol}} the user has to choose first between the
basis help and the help provided by \code{spam}.\cr

Sparseness is used when handling large matrices. Hence, care has been
used to provide efficient and fast routines. Essentially, the functions
do not transform the sparse structure into full matrices to use standard
(available) functionality, followed by a back transform. We agree, more
operators, functions, etc. should eventually be implemented.


The packages \code{fields} and \code{spdep} use \code{spam} as a
required package.  \cr
  
}
\references{\url{www.math.uzh.ch/furrer/software/spam/}
}
\seealso{See \code{\link{spam.class}} for a detailed class description,
  \code{\link{spam}} and \code{\link{spam.ops}} for creation,
  coercion and algebraic operations.\cr%
%  
\code{demo(package='spam')} lists available demos.\cr% 
Related packages are \code{\link[fields]{fields}},
  \code{\link[Matrix]{Matrix}} and
  \code{\link[SparseM]{SparseM.ontology}}.
}
\examples{
\dontrun{
## History of changes
file.show(system.file("NEWS", package = "spam"))
}
}
\author{Reinhard Furrer}
\keyword{documentation}
\keyword{package}
