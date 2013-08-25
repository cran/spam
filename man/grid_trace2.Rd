% This is file ../spam/man/grid_trace2.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{grid_trace2}
\alias{grid_trace2}
\title{Two trace plots and a scatter plot.}
\description{For two (MCMC) chains of the same length trace plots
  and a scatter plot are drawn. }
\usage{
grid_trace2(chain1, chain2 = NULL,
                   chain1_lab = '',
                   chain2_lab = '',
                   chain1_yaxis_at = NULL,
                   chain2_yaxis_at = NULL,
                   log = FALSE,
                   draw = TRUE)
}
\arguments{
\item{chain1}{A numeric vector or alternatively a matrix with two columns.}
\item{chain2}{A numeric vector of same length as \code{chain1}. (Only
  needed if \code{chain1} is specified as vector).}
\item{chain1_lab}{Label of chain 1 (character of length 1). }
\item{chain2_lab}{Label of chain 2 (character of length 1). }
\item{chain1_yaxis_at}{Points at which tick-marks are drawn for chain 1.}
\item{chain2_yaxis_at}{Points at which tick-marks are drawn for chain 2.}
\item{log}{Logical, should the date be log transformed?}
\item{draw}{Logical, should the returned grob object be drawn?}
}
\value{A grob object.}
%\details{}
%\references{}
\seealso{grid_zoom}
\examples{
grid_trace2(runif(500),runif(500),
            chain1_yaxis_at = seq(.2, 1, by = .2),
            chain1_lab = 'chain1', chain2_lab = 'chain2')
}
\author{Florian Gerber}
\keyword{hplot}
\keyword{grid}
