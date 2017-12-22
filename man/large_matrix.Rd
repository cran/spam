% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/diag.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{large_matrix}
\alias{large_matrix}
\alias{large matrix}
\alias{large_matries}
\alias{large matries}
\alias{64bit}
\alias{64bits}
\alias{long vector}
\alias{long vectors}
\alias{spam64}

\title{Large 64-bit matrices require the R package \pkg{spam64}}
\description{The R package \pkg{spam} can handle sparse matrices with up
  to 2^31-1 non-zero elements. For matrices with more non-zero elements
  it is necessary to load the \pkg{spam64} package in addition. 
}
\details{With the help of the R package \pkg{dotCall64} spam interfaces
  either the compiled code with 32-bit integers provided in
  \pkg{spam} or the compiled code with 64-bit integers provided in
  \pkg{spam64}.
}
\references{
F. Gerber, K. Moesinger, R. Furrer (2017),
Extending R packages to support 64-bit compiled code: An illustration with spam64 and GIMMS NDVI3g data,
Computer & Geoscience 104, 109-119, https://doi.org/10.1016/j.cageo.2016.11.015.
}
\seealso{
 \code{\link[spam64]{spam64-package}},  \code{\link[dotCall64]{dotCall64}}. 
}
\examples{
\dontrun{
## the following matrices are very large, and hence, 
## require much memory and cpu time.
library("spam64")
s1 <- spam(1, ncol=2^30)        # 32-bit matrix
s1

s2 <- cbind(s1, s1)             # 64-bit matrix
s2

s3 <- spam(1, ncol=2^31)        # 64-bit matrix
s3
}
}
\author{Reinhard Furrer, Florian Gerber, Kaspar Moesinger, Daniel Gerber}
\keyword{array}
\keyword{algebra}
