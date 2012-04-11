# This is file ../spam0.29-0/R/helper.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# written and maintained by Reinhard Furrer.
     



########################################################################
########################################################################
# a few nice helper functions:


bandwidth <- function(A) {
  if (!is.spam(A)) {
    warning("Matrix not 'spam' object. Coerced to one")
    A <- as.spam(A)
  }
  ret <- .Fortran("getbwd",A@dimension[1],A@entries,A@colindices,
                  A@rowpointers,low=integer(1),upp=integer(1),
                 NAOK = !.Spam$safemode[3], DUP=FALSE, PACKAGE = "spam")
  return(c(ret$low,ret$upp))
}
                  
  

bdiag.spam <- function(...){
  nargs <- nargs()
  if (nargs == 0)     return( NULL)
  args <- list(...)
  args[which(sapply(args, is.null))] <- NULL

  if (nargs == 1)     return( args[[1]])
  if (nargs == 2) {
    # Classical case, concatenate two matrices
    A <- args[[1]]
    B <- args[[2]]
    if(!is.spam(A))
      A <- as.spam(A)
    if(!is.spam(B))
      B <- as.spam(B)
    dimA <- A@dimension
    lenA <- length(A@entries)

    A@entries <- c(A@entries,B@entries)
    A@colindices <- c(A@colindices,B@colindices+dimA[2])
    A@rowpointers <- c(A@rowpointers,B@rowpointers[-1]+lenA) 
    A@dimension <-  dimA+B@dimension
    return(A)
  } else {
    # "recursive" approach only, e.g. no checking
    tmp <- bdiag.spam( args[[1]], args[[2]])
    for ( i in 3:nargs)
      tmp <- bdiag.spam( tmp, args[[i]])
    return( tmp)
  }
}


adjacency.landkreis <- function(loc)
  # this reads the germany graph file provide by
  # loc <- "http://www.math.ntnu.no/~hrue/GMRF-book/germany.graph"
  # or
  # loc <- system.file("demodata/germany.graph", package="INLA")
  # 
  {
    n <- as.numeric( readLines(loc, n=1))

    nnodes <- nodes <- numeric( n)

    adj <- list()
    for (i in 1:n) {
      tmp <- as.numeric(scan(loc, skip=i, nlines=1, quiet=T, what=list(rep("",13)))[[1]])
      nodes[i] <- tmp[1]
      nnodes[i] <- tmp[2]
      adj[[i]] <- tmp[-c(1:2)]
    }
    
    adj <- adj[ order(nodes)]
    nnodes <- nnodes[ order(nodes)]
    A <- spam(0)
    A@colindices <- as.integer( unlist(adj)+1)
    A@rowpointers <- as.integer( c(1,cumsum(lapply(adj, length))+1))
    A@entries <- rep(1, length(unlist(adj)))
    A@dimension <- as.integer( c(n, n))
    return(A)
  }

map.landkreis <- function(data, col=NULL, zlim=range(data), add=FALSE, legendpos=c( 0.88,0.9,0.05,0.4))
# This is a stripped-down version of the function provided by the INLA package.
# Added color argument, changed 'append' to 'add'.
# Legend is tuned for a mai=rep(0,4) call 
{
  npoly <- length(germany)
  ymax <- ymin <- xmax <- xmin <- 1:npoly

  if (length(data)!=npoly)
    stop('data has wrong length')
  
  if (is.null(col)) {
    if (exists('tim.colors'))
      col <- tim.colors(64)
    else
      col <- gray(seq(.05,to=0.95,length=64))
  }
  ncol <- length(col)
  polycol <- col[round(((data-zlim[1])/diff(zlim)+1e-6)*(ncol-1))+1]
  
  for(i in 1:length(germany)) {
    xmin[i] <- min(germany[[i]][,2],na.rm=T)
    xmax[i] <- max(germany[[i]][,2],na.rm=T)
    ymin[i] <- min(germany[[i]][,3],na.rm=T)
    ymax[i] <- max(germany[[i]][,3],na.rm=T)
  }


  if (!add)
    plot(c(min(xmin),max(xmax)),c(min(ymin),max(ymax)), type="n", axes=F, xlab="", ylab="")
  for(k in npoly:1)
    polygon(germany[[k]][,2],germany[[k]][,3],col=polycol[k])
  if (exists('image.plot'))
    image.plot(as.matrix(data), zlim=zlim, legend.only=T, smallplot=legendpos, cex=.2, col=col)

  invisible()
}
