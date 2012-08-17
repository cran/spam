# This is file ../spam0.29-2/R/rowcolstats.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# written and maintained by Reinhard Furrer.
     



rowSums.spam <- function(x,...) {
  return( .Fortran("rowsums",
                   dcheck(x@entries), icheck(x@colindices), icheck(x@rowpointers),
                   x@dimension[1],
                   rs=vector("double",x@dimension[1]),
                   NAOK=!.Spam$safemode[3], DUP = FALSE, PACKAGE="spam")$rs)
  
}

colSums.spam <- function(x,...) {
  return( .Fortran("colsums",
                   dcheck(x@entries), icheck(x@colindices), icheck(x@rowpointers),
                   x@dimension[1],
                   cs=vector("double",x@dimension[2]),
                   NAOK=!.Spam$safemode[3], DUP = FALSE, PACKAGE="spam")$cs)
}

rowMeans.spam <- function(x,...) {
  return( .Fortran("rowmeans",
                   dcheck(x@entries), icheck(x@colindices), icheck(x@rowpointers),
                   x@dimension[1],x@dimension[2],
                   as.logical(.Spam$structurebased),
                   rm=vector("double",x@dimension[1]),
                   NAOK=!.Spam$safemode[3], DUP = FALSE, PACKAGE="spam")$rm)
}

colMeans.spam <- function(x,...) {
   return( .Fortran("colmeans",
                dcheck(x@entries), icheck(x@colindices), icheck(x@rowpointers),
                x@dimension[1],x@dimension[2],
                as.logical(.Spam$structurebased),
                cm=vector("double",x@dimension[2]),vector("integer",x@dimension[2]),
                NAOK=!.Spam$safemode[3], DUP = FALSE, PACKAGE="spam")$cm)
}



setMethod("rowSums","spam",rowSums.spam)
setMethod("colSums","spam",colSums.spam)
setMethod("rowMeans","spam",rowMeans.spam)
setMethod("colMeans","spam",colMeans.spam)
