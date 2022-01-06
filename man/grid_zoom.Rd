% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/grid_zoom.Rd.                       %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







\name{grid_zoom}
\alias{grid_zoom}
\title{grid_zoom}
\description{This function takes a grob object (e.g. created with
  package grid) and adds a zoom window. }
\usage{
grid_zoom(inputGrob = pointsGrob(runif(200),runif(200)),
          inputViewport = viewport(name='main'),
          x = 'topleft', y, just,
          ratio = c(.3,.4), zoom_xlim, zoom_ylim,
          rect = TRUE, rect_lwd = 1, rect_fill = 'gray92',
          draw =TRUE, zoom_fill = 'white',
          zoom_frame_gp = gpar(lwd = 1),
          zoom_gp = NULL, zoom_xaxis = xaxisGrob(main = FALSE),
          zoom_yaxis = NULL)
}
\arguments{
  \item{inputGrob}{A grob object, e.g created with package grid.}
  \item{inputViewport}{Viewport related to \code{inputGrob}.}
  \item{x}{Specifies the \code{x} coordinate of the zoom
    window. Alternatively it can be set to 'topleft', 'topright',
    'bootmleft' or 'bootmright'}
  \item{y}{Specifies the \code{y} coordinate of the zoom
    window. }
  \item{just}{Specifies the justification of the zoom window. }
  \item{ratio}{Specifies size of the zoom window relative to the main window. }
  \item{zoom_xlim}{Specifies xlim value of the zoom window. }
  \item{zoom_ylim}{Specifies ylim value of the zoom window. }
  \item{rect}{Logical, if TRUE a rectangle of the zoom region is draw in the main window. }
  \item{rect_lwd}{lwd of the rectangle. }
  \item{rect_fill}{fill of the rectangle. }
  \item{draw}{logical, if TRUE the returned grob object is also drawn.}
  \item{zoom_fill}{fill color of the zoom window.}
  \item{zoom_frame_gp}{gpar() of the frame of the zoom window.}
  \item{zoom_gp}{gpar() of the inputGrob in the zoom viewport.}
  \item{zoom_xaxis}{xaxisGrob() to draw for the zoom window.}
  \item{zoom_yaxis}{yaxisGrob() to draw for the zoom window.}
}
\value{A grob object.}
\details{
A zoom plot does only make sense if all objects of the \code{inputGrob}
are specified in \code{native} units. Additional caution me be require
for certain grobs: e.g. a zoom of a circleGrob() is problematic if the x
and y axis are stretched by a different amount.

}
%\references{}
\seealso{grid_trace2}
\examples{
require(grid)
## -- Example 1 --
set.seed(133)
grid_zoom(inputGrob = pointsGrob(runif(200), runif(200)),
          inputViewport = viewport(name = 'main'),
          zoom_xlim = c(.2, .3), zoom_ylim = c(.2, .3))


## -- Example 2 --
## initial plot
grid.newpage()
vp <- viewport(width=.8, height=.8, clip='on')
gt <- gTree(children=gList(polylineGrob(x=c((0:4)/10, rep(.5, 5), (10:6)/10, rep(.5, 5)),
              y=c(rep(.5, 5), (10:6/10), rep(.5, 5), (0:4)/10),
              id=rep(1:5, 4), default.units='native',
              gp=gpar(col=1:5, lwd=3)),
              pointsGrob(runif(1000), runif(1000),pch='.', gp=gpar(cex=3)),
              rectGrob(gp=gpar(lwd=3))))
pushViewport(vp)
grid.draw(gt)

## plot with zoom window
grid.newpage()
grid_zoom(inputGrob = gt,
          inputViewport = vp,
          x='topright', zoom_xlim=c(.6,.73), zoom_ylim=c(.3,.43),ratio=.4,
          zoom_xaxis = NULL, zoom_gp = gpar(cex=3))
}
\author{Florian Gerber <florian.gerber@math.uzh.ch>}
\keyword{zoom_grid}
\keyword{zoom}
\keyword{plot}
\keyword{grid}
