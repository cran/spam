# This is file ../spam0.28-0/R/plotting.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# written and maintained by Reinhard Furrer.
     










image.spam <- 
function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)),
    z, zlim = range(z), xlim = range(x), ylim = range(y),
    col = heat.colors(12), add = FALSE, xaxs = "i", yaxs = "i",
    xlab, ylab, breaks, oldstyle = FALSE,cex=NULL, ...)
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                if (is.null(dim(x)))
                  stop("argument must be matrix-like")
                z <- x
                x <- seq(0, 1, len = nrow(z))
            }
            if (missing(xlab))
                xlab <- ""
            if (missing(ylab))
                ylab <- ""
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        xn <- deparse(substitute(x))
        if (missing(xlab))
            xlab <- paste(xn, "x", sep = "$")
        if (missing(ylab))
            ylab <- paste(xn, "y", sep = "$")
        y <- x$y
        x <- x$x
    }
    else {
        if (missing(xlab))
            xlab <- if (missing(x))
                ""
            else deparse(substitute(x))
        if (missing(ylab))
            ylab <- if (missing(y))
                ""
            else deparse(substitute(y))
    }
    spamversion <- (prod(z@dimension) > .Spam$imagesize)

    if (any(!is.finite(x)) || any(!is.finite(y)))
      stop("'x' and 'y' values must be finite and non-missing")
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    if (!is.spam(z))    stop("'z' must be a matrix")
    if (spamversion) {
      xx <- x
      yy <- y
    }
      if (length(x) > 1 && length(x) == nrow(z)) {
        dx <- 0.5 * diff(x)
        x <- c(x[1] - dx[1], x[-length(x)] + dx, x[length(x)] +
               dx[length(x) - 1])
      }
      if (length(y) > 1 && length(y) == ncol(z)) {
        dy <- 0.5 * diff(y)
        y <- c(y[1] - dy[1], y[-length(y)] + dy, y[length(y)] +
               dy[length(y) - 1])
      }
    
    if (!spamversion) {
      zvals <- as.matrix(z)
      zvals[zvals==0] <- NA
     } else zvals <- z@entries
    if (missing(breaks)) {
        nc <- length(col)
        if (!missing(zlim) && (any(!is.finite(zlim)) || diff(zlim) < 0))
            stop("invalid z limits")
        if (diff(zlim) == 0)
            zlim <- if (zlim[1] == 0) {
                c(-1, 1)
            } else zlim[1] + c(-0.4, 0.4) * abs(zlim[1])
        zvals <- (zvals - zlim[1])/diff(zlim)
        zi <- if (oldstyle) {
            floor((nc - 1) * zvals + 0.5)
        } else floor((nc - 1e-05) * zvals + 1e-07)
        zi[zi < 0 | zi >= nc] <- NA
    }
    else {
        if (length(breaks) != length(col) + 1)
            stop("must have one more break than colour")
        if (any(!is.finite(breaks)))
            stop("breaks must all be finite")
        zi <- .C("bincode", as.double(zvals), length(zvals), as.double(breaks),
            length(breaks), code = vector("integer",length(zvals)), (TRUE),
            (TRUE), nok = TRUE, NAOK = TRUE, DUP = FALSE, PACKAGE = "base")$code -
            1
    }
    if (!add)
        plot(NA, NA, xlim = xlim, ylim = ylim, type = "n", xaxs = xaxs,
            yaxs = yaxs, xlab = xlab, ylab = ylab, ...)
    if (spamversion) {
      if (length(xx) != nrow(z) || length(yy) != ncol(z))
        stop("dimensions of z are not length(x) times length(y)")
    }else{
      if (length(x) <= 1)
        x <- par("usr")[1:2]
      if (length(y) <= 1)
        y <- par("usr")[3:4]
      if (length(x) != nrow(z) + 1 || length(y) != ncol(z) + 1)
        stop("dimensions of z are not length(x)(+1) times length(y)(+1)")
    }
# for small matrices, we transform them into regular ones.
    if (!spamversion) {
      .Internal(image(as.double(x), as.double(y), as.integer(as.matrix(zi)),col))
    } else {
      if (missing(cex)) {
        warning("default value for 'cex' in 'image' might be a bad choice", call.=FALSE)
        cex <- 1
      }
      points( xx[rep.int((1:nrow(z)),diff(z@rowpointers))], yy[z@colindices],
             pch='.', cex=cex*.Spam$cex/sum(z@dimension),
             col=col[zi+1])
    }
    box()
  }

display.spam <- function(x,col=c("gray","white"),xlab="column",ylab="row", cex=NULL,
                       main="",...)
{
  nrow <- x@dimension[1]
  ncol <- x@dimension[2]
  
# For small matrices, we transform them into regular ones and use the image.default
# routine.  
  if (prod(nrow,ncol) < .Spam$imagesize) {
    z <- vector("double", prod(nrow,ncol))
    dim(z) <- c(nrow,ncol)
    z[cbind(rep.int(nrow:1,diff(x@rowpointers)),x@colindices)] <- -1
    image.default(x=1:ncol,y=-(nrow:1),t(z),
                  axes=FALSE, col=col, xlab=xlab, ylab=ylab, ...) 
  } else {
    if (missing(cex)) {
      warning("default value for 'cex' in 'display' might not be the optimal choice", call.=FALSE)
      cex <- 1
    }
    plot( x@colindices, rep.int(-(1:nrow),diff(x@rowpointers)), pch='.', cex=cex*.Spam$cex/(ncol+nrow),
         col=col[1],xlab=xlab,ylab=ylab,axes=FALSE,
         ylim=c(-nrow,-0)-.5,xlim=c(0,ncol)+.5,xaxs = "i", yaxs = "i",...)
  }
  # Adjust axes labels.
  axis(1,pretty(1:ncol), ...)
  axis(2,pretty(-(nrow:1)),labels=rev(pretty(1:nrow)), ...)
  box()
}



plot.spam <- function(x,y,xlab=NULL,ylab=NULL,...)
{
  lab <- deparse(substitute(x))
  #only a few cases are considered
  # 1st case, a colum vector only
  if (ncol(x)==1) {
    x <- c(x)
    return( plot(x,...))
  }
  # 2nd case a matrix
  tmp <- x[,1:2, drop=FALSE] # extract the first two columns
  plot(c( tmp[,1]), c(tmp[,2]),
       xlab=ifelse(missing(xlab),paste(lab,'[,1]',sep=''),xlab),
       ylab=ifelse(missing(ylab),paste(lab,'[,2]',sep=''),ylab),...)
}

setGeneric("image", function(x, ...) standardGeneric("image")) 
setMethod("image","spam",function(x,cex=NULL,...){image.spam(x,cex=cex,...)})

setGeneric("display",function(x,...)standardGeneric("display"))
setMethod("display","spam",display.spam)

setMethod("plot", signature(x="spam",y="missing"), plot.spam)
setMethod("plot", signature(x="spam",y="spam"),
          function(x,y,...) {
            warning("'plot' with two 'spam' objects is not implemented",call.=FALSE)
            })


