% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/grid_trace2.Rd.                     %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

\name{grid_trace2}
\alias{grid_trace2}
\title{Two trace plots and a scatter plot.}
\description{For two (MCMC) chains of the same length trace plots
  and a scatter plot are drawn. }
\usage{
grid_trace2(chain1, chain2 = NULL,
            xlim = NULL, ylim1 = NULL, ylim2=NULL,
            chain1_lab = "", chain2_lab = "", main = "",
            chain1_yaxis_at = NULL, chain2_yaxis_at = NULL,
            log = FALSE,
            cex_points = unit(0.2, "mm"),
            cex_main = unit(1.2, "mm"),
            lwd_lines = unit(0.1, "mm"),
            lwd_frame = unit(0.8, "mm"),
            draw = TRUE)
}
\arguments{
\item{chain1}{Numeric vector or a matrix with two columns.}
\item{chain2}{Numeric vector of same length as \code{chain1}. (Only used if \code{chain1} is specified as vector).}
\item{xlim}{x axis limits of both chains (numeric vector of length 2). }
\item{ylim1}{y axis limits of chain 1 (numeric vector of length 2).}
\item{ylim2}{y axis limits of chain 2 (numeric vector of length 2).}
\item{chain1_lab}{Label of chain 1 (character of length 1). }
\item{chain2_lab}{Label of chain 2 (character of length 1). }
\item{main}{Title (character of length 1). }
\item{chain1_yaxis_at}{Points at which tick-marks are drawn for chain 1.}
\item{chain2_yaxis_at}{Points at which tick-marks are drawn for chain 2. }
\item{log}{Logical, should the date be log transformed?}
\item{cex_points}{Size of points in scatter plot. }
\item{cex_main}{Size of the title font. }
\item{lwd_lines}{Line width of trace plots. }
\item{lwd_frame}{Line width of frames. }
\item{draw}{Logical, should the returned grob object be drawn?}
}
\value{A grob object.}
\details{The figure is optimized for a plot of the format
  x11(width = 10, height = 3).}
%\references{}
\seealso{grid_zoom}
\examples{
grid_trace2(runif(500), runif(500),
            chain1_yaxis_at = seq(.2, 1, by = .2),
            chain1_lab = "chain1", chain2_lab = "chain2",
            main = "Uniform",
            lwd_lines = unit(.5, "mm"))
}
\author{Florian Gerber <florian.gerber@math.uzh.ch>}
\keyword{plot}
\keyword{trace plot}
\keyword{grid}
