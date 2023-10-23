% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/spam-package.Rd.                    %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\name{spam-package}
\alias{spam-package}
\alias{SPAM}
\alias{Spam}
\alias{. SPAM .}
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
function. For example, \code{help(chol.spam)} calls the help
directly. We deliberately avoided aliases according to analogue helps from the
base package.\cr

Sparseness is used when handling large matrices. Hence, care has been
used to provide efficient and fast routines. Essentially, the functions
do not transform the sparse structure into full matrices to use standard
(available) functionality, followed by a back transform. We agree, more
operators, functions, etc. should eventually be implemented.


The packages \code{fields} and \code{spam} are closely linked.  \cr

}
\references{Reinhard Furrer, Stephan R. Sain (2010).
  "spam: A Sparse Matrix R Package with Emphasis on MCMC
  Methods for Gaussian Markov Random Fields.",
        \emph{Journal of Statistical Software}, 36(10), 1-25,
        \doi{10.18637/jss.v036.i10}.\cr
Florian Gerber, Reinhard Furrer (2015).
        "Pitfalls in the Implementation of Bayesian Hierarchical Modeling of Areal Count Data: An Illustration Using BYM and Leroux Models.",
        \emph{Journal of Statistical Software}, Code Snippets, 63(1), 1-32,
        \doi{10.18637/jss.v063.c01}.\cr
F. Gerber, K. Moesinger, R. Furrer (2017),
        "Extending R packages to support 64-bit compiled code: An illustration with spam64 and GIMMS NDVI3g data.",
        \emph{Computer & Geoscience} 104, 109-119, \doi{10.1016/j.cageo.2016.11.015}."
}

\seealso{See \code{\link{spam.class}} for a detailed class description,
  \code{\link{spam}} and  \code{\link{spam.ops}} for creation,
  coercion and algebraic operations. \code{\link{options}}.%\cr%
%
%\code{demo(package='spam')} lists available demos.\cr%
%Related packages are \code{\link[fields]{fields}},
%  \code{\link[Matrix]{Matrix}} and
%  \code{SparseM.ontology}.
}
\examples{
## Citations:
citation('spam')
citation('spam', auto=TRUE)


## History of changes
\dontrun{
file.show(system.file("NEWS.md", package = "spam"))
}
}
\author{Reinhard Furrer, with the help of Florian Gerber, Kaspar
  Moesinger and many others. \cr
  Some Fortran routines were written by Youcef Saad, Esmond G. Ng, Barry W. Peyton, Joseph W.H. Liu, Alan D. George.
}
\keyword{documentation}
\keyword{package}
