% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/Oral.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   




\name{Oral}
\alias{Oral}
\alias{Oral.spam}
\alias{oral.spam}
\docType{data}
\title{Oral Cavity Cancer}
\description{Oral cavity cancer counts in 544 districts in Germany over 1986-1990.}
\format{\code{Oral} is a dataframe with 3 columns. 
\describe{
  \item{Y}{observed counts}
  \item{E}{expected counts}
  \item{SMR}{standardized mortality ratios}}
\code{germany} is a
list of 544 elements, each describing an individual polygon of the district.
}
\details{The expected counts depend on the number of people in the
  region and their age distribution.\cr
The regions are ordered according the supplied polygon description and
adjacency graph.\cr

There is a similar dataset \code{data(Germany)}  with larynx cancer cases from the package
  \pkg{INLA} as well, with an additional smoking covariate.
}
\source{The data is available from the package \pkg{INLA}
  distributed from \url{www.r-inla.org} or from\cr
  \url{http://www.math.ntnu.no/~hrue/GMRF-book/oral.txt}
}
\references{
Knorr-Held, L. and Rasser, G. (2000) Bayesian Detection of Clusters and Discontinuities in Disease Maps,
     \emph{Biometrics},
     56,
     13--21.
}
\seealso{\code{\link{germany.plot}}.}
\keyword{datasets}
