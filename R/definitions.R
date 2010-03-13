# This is file ../spam0.21-0/R/definitions.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# written and maintained by Reinhard Furrer.
     



todo <- function() help( "spam.todo")
spam.history <- function() help("spam.history")


dcheck <- function(x) if (.Spam$safemode[1]) as.double(x) else x
icheck <- function(x) if (.Spam$safemode[1]) as.integer(x) else x

validspamobject <- function(object) {

  if (.Spam$safemode[2]){ 
    if(!identical(length(object@dimension), int2) ){
      print(object@dimension)
      return("invalid dimension attribute")
    }
    else{
      nrow <- object@dimension[1]
      ncol <- object@dimension[2]
    }
    if (any(!is.double(object@entries)))
      return("matrix entries need to be of double mode")
    if (any(!is.finite(object@entries)))
      return("'NA/NaN/Inf' not allowed")
    if(!identical(length(object@entries),length(object@colindices)))
      return("entries and column indices don't have equal lengths")
    if(any(object@colindices < 1) || any(object@colindices > ncol)) 
      return("column indices exceeds dimension bounds")
    if(any(object@rowpointers < 1))
      return("some elements of row pointes are <= 0")
    if(any(diff(object@rowpointers)<0))
      return("row pointers are not monotone increasing")
    diffcolindices <- diff(object@colindices)     # positive values within each row
    if (all(diff(object@rowpointers)>1) && length(diffcolindices)>0)   # only if we have multiple values
      if (identical( nrow, int1)) {
        if   ( any(diffcolindices<1))
          return("column indices are not ordered")
      } else {
        if ( any(diffcolindices[-(object@rowpointers[2:nrow]-1)]<1))
          return("column indices are not ordered")
      }
    if(object@rowpointers[length(object@rowpointers)] != length(object@entries)+1)
      return("last element of row pointers doesn't conform")
    if(length(object@rowpointers) != nrow+1)
      return("row pointers has wrong number of elements")
    if(length(object@entries) < 1 || length(object@entries) > prod(object@dimension))
      return("too few or too many entries")
  }
  TRUE
}

setClass("spam",representation(entries="numeric",
                               colindices="integer",
                               rowpointers="integer",
                               dimension="integer"),
         prototype = prototype(entries=as.double( 0),
                               colindices=as.integer(1),
                               rowpointers=as.integer(c(1,2)),
                               dimension=as.integer(c(1,1))),
         validity  = validspamobject)

setMethod("initialize", "spam", 
          function(.Object,
                   entries    = 0,         # by default a 1x1 zero matrix.
                   colindices = as.integer( rep(1, length( entries ))),  # or a nx1 matrix
                   rowpointers= as.integer( 1:(length( entries )+1)),   # with n=length(ra)
                   dimension  = as.integer( c(length( rowpointers )-1,max( 1,colindices ))))
          {
            # a specific "degenerate" case:
            if (identical(length(entries),int0)) {  # e.g., induced by rep(1,0)
              warning("While initializing, empty 'spam' object coerced to zero 'spam' matrix",call.=FALSE)
              entries <- 0
              colindices <- int1
              rowpointers <- c(int1,int2)
              dimension <- c(int1,int1)
            }
#            if (rowpointers[ length(rowpointers)] ==1) {  # e.g., zero matrix
#              rowpointers[-1] <- as.integer(2)
#              colindices <- as.integer(1)
#            }
            .Object@entries     <- entries
            .Object@colindices  <- colindices
            .Object@rowpointers <- rowpointers
            .Object@dimension   <- dimension
            validObject(.Object)
            .Object
          })


print.spam <- function(x,...) {
  if (prod(x@dimension) < .Spam$printsize) {
    print(as.matrix.spam(x),...)
  } else {
    if ( (length(x@entries)==1)  & (x@entries[1]==0)) {
      cat('Zero matrix of dimension ',x@dimension[1],'x',
                x@dimension[2],'.\n',sep='', fill=TRUE)
    }else {
      cat('Matrix of dimension ',x@dimension[1],'x',
                x@dimension[2],' with (row-wise) nonzero elements:\n',sep='', fill=TRUE)
      print(x@entries,...)
    }
  }
  cat("Class 'spam'\n")
  invisible(NULL)
}

summary.spam <- function(object,...) {
            nz <- length(object@entries)
            dens <- nz/prod(object@dimension)*100
            cat("Matrix object of class 'spam' of dimension ",object@dimension[1],'x',
                object@dimension[2],',\n',sep='')
            cat('    with ',nz,' (row-wise) nonzero elements.\n',sep='')
            cat('    Density of the matrix is ',signif(dens,3),'%.\n',sep='')
            cat("Class 'spam'\n")
            invisible(list(nnz=nz, density=dens))
          }

"dim<-.spam" <- function(x,value) {
  if ( (min(value)<1 ) || any(!is.finite(value)))
    stop("dims should be postive integers.")
  if (!identical( length(value), int2)) stop("dims should be of length 2.")
  dimx <- x@dimension
  last <- value[1]+1

  # In three steps:
  #  1) Address col truncation
            # to safe time, we also take into account if we have fewer or equal rows
  #  2) Augment rows
  #  3) if fewer rows and more columns, truncate
  # In any case, dimensions are fixed at the end.
  
  # If fewer cols required, we run reducedim
  if (dimx[2]>value[2]){
#     subroutine reducedim(a,ja,ia,eps,bnrow,bncol,k,b,jb,ib)
    z <- .Fortran("reducedim",
                  oldra=dcheck(x@entries),
                  oldja=x@colindices,
                  oldia=x@rowpointers,
                  eps=.Spam$eps,
                  as.integer(min(value[1],dimx[1])),as.integer(value[2]),
                  nz=as.integer(1),
                  entries=vector("double",length(x@entries)),
                  colindices=vector("integer",length(x@entries)),
                  rowpointers=vector("integer",last),
                  NAOK = !.Spam$safemode[3], DUP=FALSE, PACKAGE = "spam")
    if (identical(z$nz,int1) )
      return(new("spam",rowpointers=c(int1,rep.int(int2,as.integer(value[1]))),
                 dimension=as.integer(value)))

    nz <- z$nz-1
    slot(x,"entries",check=FALSE) <- z$entries[1:nz]
    slot(x,"colindices",check=FALSE) <- z$colindices[1:nz]
    slot(x,"rowpointers",check=FALSE) <- z$rowpointers[1:min(last,dimx[1]+1)]
  }
  # augment rows
  if  (dimx[1]<value[1]){
    slot(x,"rowpointers",check=FALSE) <- c(x@rowpointers,rep.int(
                           x@rowpointers[length(x@rowpointers)],value[1]-dimx[1]))
  }
  # special case: fewer rows and more columns, truncate
  if ((dimx[1]>value[1])&(dimx[2]<value[2])) {
    lastelement <- (x@rowpointers[last]-1)
    slot(x,"entries",check=FALSE) <- x@entries[1:lastelement]
    slot(x,"colindices",check=FALSE) <- x@colindices[1:lastelement]
    slot(x,"rowpointers",check=FALSE) <- x@rowpointers[1:last]
  }
        
  slot(x,"dimension",check=FALSE) <- as.integer(value)                  
  return(x)

}

setMethod("show","spam",  function(object) {
  if (prod(object@dimension) < .Spam$printsize) {
    print(as.matrix.spam(object))
  } else {
    if ( identical(length(object@entries),int1)  & identical(object@entries[1],int0)) {
      cat('Zero matrix of dimension ',object@dimension[1],'x',
                object@dimension[2],'.\n',sep='')
    }else {
      cat('Matrix of dimension ',object@dimension[1],'x',
                object@dimension[2],' with (row-wise) nonzero elements:\n',sep='')
      print(object@entries)
    }
  }
  cat("Class 'spam'\n")
  invisible(object)
})

setMethod("print","spam", print.spam)
setMethod("summary","spam",summary.spam)
setMethod("length","spam",function(x) x@rowpointers[x@dimension[1]+1]-1 ) # equivalent to length(x@entries) 
setMethod("dim",   "spam",function(x) x@dimension )

setMethod("length<-","spam",function(x,value) stop("operation not allowed on 'spam' object") )
setMethod("dim<-",   "spam", get("dim<-.spam"))


setMethod("c","spam", function(x,...,recursive=TRUE){
  dimx <- x@dimension
  cx <- .Fortran("spamcsrdns",
                 nrow=dimx[1],
                 entries=dcheck(x@entries),
                 colindices=x@colindices,
                 rowpointers=x@rowpointers,
                 res=vector("double",prod(dimx)),  
                 NAOK=!.Spam$safemode[3],DUP=FALSE,PACKAGE = "spam")$res
  if (length( list(...)) < 1)
    return( cx)
  else
    c( cx,c(...,recursive),recursive)
})

########################################################################
# diag and derivatives
"diag.spam" <- function(x=1, nrow, ncol)  {
  if (is.spam(x)) return( diag.of.spam( x, nrow, ncol))

  if (is.array(x) && length(dim(x)) != 1)
    stop("first argument is array, but not matrix.")

  if (missing(x))
    n <- as.integer(nrow)
  else if (length(x) == 1 && missing(nrow) && missing(ncol)) {
    n <- as.integer(x)
    x <- 1
  }
  else n <- length(x)
  if (!missing(nrow))
    n <- as.integer(nrow)
  if(missing(ncol))
    ncol <- n
  p <- as.integer(ncol)

  m <- min(n, p)

  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- vector("double", m)
  newx@entries[1:m] <- as.double(x) 
  slot(newx,"colindices",check=FALSE) <- 1:m
  slot(newx,"rowpointers",check=FALSE) <- as.integer(c(1:m,rep(m+1,n+1-m)))
  slot(newx,"dimension",check=FALSE) <- c(n,p)
  return(newx)
}



"diag<-.spam" <-  function(x,value) {
  nrow <- x@dimension[1]
  minrc <- min( x@dimension)
  if (length(value)==1)
    value <- rep(value,minrc)
  else if (length(value)!=minrc)
    stop("replacement diagonal has wrong length")
  z <- .Fortran("setdiagmat",
                nrow = nrow,
                n = minrc,
                a = c(x@entries,double(minrc)),
                ja = c(x@colindices,integer(minrc)),
                ia = x@rowpointers,
                diag = as.double(value),
                iw = vector("integer",nrow),    # just to be sure
                info = vector("integer",nrow+1),
                NAOK = !.Spam$safemode[3],
                PACKAGE = "spam")
  nz <- z$ia[nrow+1]-1
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- z$a[1:nz]
  slot(newx,"colindices",check=FALSE) <- z$ja[1:nz]
  slot(newx,"rowpointers",check=FALSE) <- z$ia
  slot(newx,"dimension",check=FALSE) <- x@dimension
  return(newx)
}

"diag.spam<-" <- get("diag<-.spam")

"diag.of.spam" <- function(x, nrow, ncol)   {
  len <- min(x@dimension)
  return(.Fortran("getdiag",
                a = dcheck(x@entries),
                colindices =  x@colindices,
                rowpointers = x@rowpointers,
                len = len,
                diag = vector("double",len),
                NAOK=!.Spam$safemode[3],
                DUP=FALSE,
                PACKAGE = "spam"
                )$diag)
}

setMethod("diag","spam",diag.of.spam)
setMethod("diag<-","spam",get("diag<-.spam"))

########################################################################

"t.spam" <- function(x){
  dimx <- x@dimension
  nz <- x@rowpointers[dimx[1]+1]-1
  z <- .Fortran("transpose",
                n=dimx[1],m=dimx[2],
                a=dcheck(x@entries),ja=x@colindices,ia=x@rowpointers,
                entries=vector("double",nz),colindices=vector("integer",nz),
                rowpointers=vector("integer",dimx[2]+1),
                NAOK=!.Spam$safemode[3],
                DUP=FALSE,
                PACKAGE = "spam")
  t.x <- new("spam")
  slot(t.x,"entries",check=FALSE) <- z$entries[1:nz]
  slot(t.x,"colindices",check=FALSE) <- z$colindices[1:nz]
  slot(t.x,"rowpointers",check=FALSE) <- z$rowpointers
  slot(t.x,"dimension",check=FALSE) <- dimx[2:1]
  return( t.x)
}
setMethod("t","spam",t.spam)

########################################################################

"is.spam" <- function(x) is(x,"spam")

"as.spam" <- function(x,  eps = .Spam$eps) stop('coercion not defined form this class')

"spam" <- function(x, nrow = 1, ncol = 1, eps = .Spam$eps) stop("argument 'x' should be of mode 'numeric' (or 'spam')")

"as.spam.spam" <- function(x, eps = .Spam$eps)  {
  if (eps<.Machine$double.eps)
    stop("'eps' should not be smaller than machine precision",call.=FALSE)
  dimx <- x@dimension
  z <- .Fortran("cleanspam",
                nrow=dimx[1],
                entries=dcheck(x@entries),
                colindices=x@colindices,
                rowpointers=x@rowpointers,
                eps=as.double(eps),
                NAOK=!.Spam$safemode[3],
                PACKAGE = "spam"
                )
  nz <- z$rowpointers[dimx[1]+1]-1
  if(nz==0) return(new("spam",rowpointers=c(int1,rep(int2,dimx[1])), dimension=dimx))

  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- z$entries[1:nz]
  slot(newx,"colindices",check=FALSE) <- z$colindices[1:nz]
  slot(newx,"rowpointers",check=FALSE) <- z$rowpointers[1:(dimx[1]+1)]
  slot(newx,"dimension",check=FALSE) <- dimx
  return(newx)
}

"as.spam.matrix" <- function(x, eps = .Spam$eps) {  
  if (eps<.Machine$double.eps)
    stop("'eps' should not be smaller than machine precision",call.=FALSE)
  dimx <- dim(x)
  nz <- length(x)
  z <- .Fortran("spamdnscsr",
                nrow=dimx[1],
                ncol=dimx[2],
                x=as.double(x),
                dimx[1],
                entries=vector("double",nz),
                colindices=vector("integer",nz),
                rowpointers=vector("integer",dimx[1]+1),
                eps=as.double(eps),
                NAOK=!.Spam$safemode[3],
                DUP = FALSE,
                PACKAGE = "spam"
                )
  nz <- z$rowpointers[dimx[1]+1]-1

  if(nz==0) return(new("spam",rowpointers=c(int1,rep(int2,dimx[1])), dimension=dimx))
                        # no nonzero values. We preserve the dimension of x
    
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- z$entries[1:nz]
  slot(newx,"colindices",check=FALSE) <- z$colindices[1:nz]
  slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
  slot(newx,"dimension",check=FALSE) <- dimx
  return(newx)
}


"as.spam.numeric" <- function(x, eps = .Spam$eps) {
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  if (any(!is.finite(x))) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[!is.finite(x)] <- 0
  }
  poselements <- abs(x)>eps
  nz <- sum(poselements)
  lx <- length(x)
  if (identical(nz,0)) # empty matrix
    return(new("spam",rowpointers=c(int1,rep.int(int2,lx)), dimension=c(lx,int1)))
    
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- as.double(x[poselements])
  slot(newx,"colindices",check=FALSE) <- rep.int(int1, nz)
  slot(newx,"rowpointers",check=FALSE) <- as.integer(cumsum(c(1, poselements)))
  slot(newx,"dimension",check=FALSE) <- c(lx,int1) 
  return(newx)
}


"as.spam.dist" <- function(x, eps = .Spam$eps) {
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  if (any(!is.finite(x))) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[!is.finite(x)] <- 0
  }
  dimx <- attr(x,"Size")
  size <- dimx*(dimx-1)/2
  z <- .Fortran("disttospam",
                nrow=dimx,
                x=as.vector(x,'double'),
                entries=vector('double',size),
                colindices=vector('integer',size),
                rowpointers=vector('integer',dimx+1),
                eps=as.double(eps),
                NAOK=!.Spam$safemode[3],
                PACKAGE = "spam"
                )
  nz <- z$rowpointers[dimx+1]-1
  if(nz==0) return(new("spam",rowpointers=c(int1,rep(int2,dimx)), dimension=c(dimx,dimx)))

  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- z$entries[1:nz]
  slot(newx,"colindices",check=FALSE) <- z$colindices[1:nz]
  slot(newx,"rowpointers",check=FALSE) <- z$rowpointers[1:(dimx+1)]
  slot(newx,"dimension",check=FALSE) <- c(dimx,dimx)
 
  return(newx)
}

"as.spam.list" <-  function(x, eps = .Spam$eps)
  spam.list(x,eps)

"spam.list" <-  function(x, nrow = 1, ncol = 1, eps = .Spam$eps) {
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  if (!is.list(x)|(length(x)<2)|(length(x)>3))
    stop("Argument 'x' needs to be a list with two or three elements")
  # two cases: list of length
  # -  two (matrix with two columns called ind* and the elements)
  # -  three (each one column called i*, j*.   

  if (identical(length(x),int2)) {
    indnr <- pmatch("ind",names(x)) 
    if (is.na(indnr)) stop("Argument 'x' needs an element called 'indices'")
    elenr <- ifelse( identical( indnr,int1), int2, int1)
    
    nz <- length( x[[elenr]])

    dimx <- dim(x[[indnr]])
    if (is.null(dimx)||(dimx[2] != 2))  stop("Indices should have two columns")
    if (dimx[1] != nz) stop("Number of indices does not match with number of elements")
    
    ir <- as.integer(x[[indnr]][,1])
    jc <- as.integer(x[[indnr]][,2])
  } else {
    inr <- pmatch("i",names(x)) 
    jnr <- pmatch("j",names(x))
    if (is.na(inr)||is.na(jnr)) stop("Argument 'x' needs elements called 'i' and 'j'")
    elenr <- c(1:3)[-c(inr,jnr)]
    nz <- length( x[[elenr]])
    
    ir <- as.integer(x[[inr]])
    jc <- as.integer(x[[jnr]])

    if ((length(ir) != nz)||(length(jc) != nz))
        stop("Number of indices does not match with number of elements")
  }
  if (identical(nz, int0))
    return(new("spam",rowpointers=c(int1,rep.int(int2,as.integer(nrow))), dimension=as.integer(c(nrow,ncol))))
  if (any( ir <= 0) || any( jc <= 0))
      stop("Indices need to be positive")
  if (any(!is.finite(x[[elenr]]))) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[[elenr]][!is.finite(x[[elenr]])] <- 0
  }
  nrow <- as.integer(ifelse(missing(nrow),max(ir),nrow))
  z <- .Fortran("triplet2csr",
                nrow=nrow, ncol=as.integer(ifelse(missing(ncol),max(jc),ncol)),
                nz,
                as.double(x[[elenr]]),ir,jc,
                entries=vector("double",nz),
                colindices=vector("integer",nz),
                rowpointers=vector("integer",nrow+1),eps,
                NAOK=TRUE, DUP = FALSE, PACKAGE = "spam"
                )
  if (identical(z$nz, int0))
    return(new("spam",rowpointers=c(int1,rep.int(int2,as.integer(nrow))), dimension=as.integer(c(nrow,ncol))))
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- z$entries
  slot(newx,"colindices",check=FALSE) <- z$colindices
  slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
  slot(newx,"dimension",check=FALSE) <- c(nrow,z$ncol)
  
  return(newx)
}


"spam.numeric" <- function(x, nrow = 1, ncol = 1, eps = .Spam$eps)  {
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  if (any(!is.finite(x))) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[!is.finite(x)] <- 0
  }
  lenx <- length( x)
  if (missing(nrow))
    nrow <- ceiling( lenx/ncol)
  else if (missing(ncol))
    ncol <- ceiling( lenx/nrow)
  dimx <- as.integer( c(nrow, ncol))
  if (lenx != prod(nrow,  ncol)) {
    if(lenx==1 && abs(x)<eps) {
      return(new("spam",rowpointers=c(int1,rep.int(int2,as.integer(nrow))), 
                 dimension=dimx))
    }
    else if(prod(nrow,ncol)%%lenx!=0)
      warning("ncol*nrow indivisable by length(x)")
      
    x <- rep(x, ceiling(prod(nrow,ncol)/lenx))
    lenx <- length( x)
  }
  z <- .Fortran("spamdnscsr",
                nrow=dimx[1],
                ncol=dimx[2],
                x=as.double(x),
                dimx[1],
                entries=vector("double",lenx),
                colindices=vector("integer",lenx),
                rowpointers=vector("integer",dimx[1]+1),
                eps=as.double(eps),
                NAOK=TRUE,
                DUP = FALSE,
                PACKAGE = "spam"
                )
  nz <- z$rowpointers[dimx[1]+1]-1

  if(nz==0) return(new("spam",rowpointers=c(int1,rep(int2,dimx[1])), dimension=dimx))
                        # no nonzero values. We preserve the dimension of x
    
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- z$entries[1:nz]
  slot(newx,"colindices",check=FALSE) <- z$colindices[1:nz]
  slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
  slot(newx,"dimension",check=FALSE) <- dimx
  return(newx)
}

# "spam.spam" <-
# function(x, nrow = 1, ncol = 1, eps = .Spam$eps)
# {
#   if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
#   x <- as.matrix.spam(x)
#   xlen <- length(x)
#   if (missing(nrow))
#     nrow <- ceiling(xlen/ncol)
#   else if (missing(ncol))
#     ncol <- ceiling(xlen/nrow)
#   if (xlen == prod(nrow, ncol))
#     dim(x) <- c(nrow, ncol)
#   else{
#     if(xlen==1 && abs(x)<eps) {
#       dimx <- c(nrow,ncol)
#       return(new("spam",entries=0,colindices=as.integer(1),
#                  rowpointers=as.integer(c(1,rep(2,dimx[1]))), 
#                  dimension=as.integer(dimx)))
#     }
#     else if(prod(nrow,ncol)%%xlen!=0)
#       warning("ncol*nrow indivisable by xlen")
#       
#     x <- rep(x, ceiling(prod(nrow,ncol)/xlen))
#     dim(x) <- c(nrow, ncol)
#   }
#   return( as.spam(x, eps=eps))
# }


if(getRversion() >= "2.10")  setOldClass(c("dist", "numeric"))


setGeneric("as.spam")
setMethod("as.spam","spam",   as.spam.spam)
setMethod("as.spam","matrix", as.spam.matrix)
setMethod("as.spam","numeric",as.spam.numeric)
setMethod("as.spam","dist",   as.spam.dist)
setMethod("as.spam","list",   { function(x,eps) spam.list(x,eps=eps)})

setGeneric("spam")
setMethod("spam","numeric",spam.numeric)
setMethod("spam","list",spam.list)
# setMethod("spam","spam",spam.spam)



triplet <- function(x, tri=FALSE){
# inverse of spam.list
  dimx <- dim(x)
  if (length(dimx)!=2) stop("'x' should be a matrix like object of dim 2")
  if (is.spam(x)) {
    return(c({ if (tri) list(i=rep(1:dimx[1],diff(x@rowpointers)),
                  j= x@colindices) else list(indices=cbind(rep(1:dimx[1],diff(x@rowpointers)),
                  x@colindices) ) }, list(values=x@entries)
             )
           )
  } else {
    return(c({ if (tri) list(i=rep.int(1:dimx[1],dimx[2]),
                             j=rep.int(1:dimx[2],rep.int(dimx[1],dimx[2]))) else
              list(indices=cbind(rep.int(1:dimx[1],dimx[2]),
                  rep.int(1:dimx[2],rep.int(dimx[1],dimx[2]))))
             } , list(values=c(x))
             )
           )
  }
}

########################################################################

if(getRversion() >= "2.5") {

    
"as.matrix.spam" <-
function(x,...){
  dimx <- x@dimension
  return(array(.Fortran("spamcsrdns",
                        nrow=dimx[1],
                        entries=dcheck(x@entries),
                        colindices=x@colindices,
                        rowpointers=x@rowpointers,
                        res=vector("double",prod(dimx)),  # numeric is double! 
                        NAOK=!.Spam$safemode[3],
                        DUP=FALSE,
                        PACKAGE = "spam"
                        )$res,
               dimx)      # we preserve dimensions
         )
  }

}else{
  
"as.matrix.spam" <-
function(x){
  dimx <- x@dimension
  return(array(.Fortran("spamcsrdns",
                        nrow=dimx[1],
                        entries=dcheck(x@entries),
                        colindices=x@colindices,
                        rowpointers=x@rowpointers,
                        res=vector("double",prod(dimx)),  # numeric is double! 
                        NAOK=!.Spam$safemode[3],
                        DUP=FALSE,
                        PACKAGE = "spam"
                        )$res,
               dimx)      # we preserve dimensions
         )
}
}


setMethod("as.matrix","spam",as.matrix.spam)

########################################################################


"complement.spam" <- function(x){
# this function returns the structure of the zeros of the spam object x. 
  nrow <- x@dimension[1]
  ncol <- x@dimension[2]
  nnz <- x@rowpointers[nrow+1]-1
  nz <- prod(nrow,ncol) - nnz
  
  # we work through special cases      
  if(nz == 0)	            return( spam(0,nrow,ncol))
  if(nnz == 1 && x@entries == 0) return( spam(1,nrow,ncol))
  # normal case, proceed to efficient function
  z <- .Fortran("notzero",
                x@colindices,
                x@rowpointers,
                nrow,
                ncol,
                as.integer(nnz),
                as.integer(nz),
                colindices = vector("integer",nz),
                rowpointers = vector("integer",nrow+1),
                NAOK=!.Spam$safemode[3],
                DUP=FALSE,
                PACKAGE="spam"
                )
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- rep.int(1.0,nz)
  slot(newx,"colindices",check=FALSE) <- z$colindices
  slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
  slot(newx,"dimension",check=FALSE) <- x@dimension 
  return(newx)
}

".spam.addsparsefull" <- function(A,B){
  # A is sparse, B is full
  if (missing(B)) return(A)
  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  pdim <- prod(nrow,ncol)
  if (is.matrix(B)) {
    if(ncol != dim(B)[2] || nrow != dim(B)[1])
      stop("non-conformable matrices")
  } else {
    if(pdim%%length(B)!=0) {
      stop("longer object length
        is not a multiple of shorter object length")
    } else  B <- rep(B,pdim %/% length(B))
  }
  if (!is.double(B[1]))  B <- as.double(B)
  return(array( .Fortran("addsparsefull",
                          nrow,dcheck(A@entries),A@colindices,
                          A@rowpointers,b=dcheck(B),PACKAGE = "spam"
                          )$b,c(nrow,ncol)))
}

".spam.subfullsparse" <- function(A,B){
  # A is sparse, B is full
  if (missing(B)) {
    A@entries <- -A@entries
    return(A)
  }
  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  pdim <- prod(nrow,ncol)
  if (is.matrix(B)) {
    if(ncol != dim(B)[2] || nrow != dim(B)[1])
      stop("non-conformable matrices")
  } else {
    if(pdim %% length(B)!=0) {
      stop("longer object length
        is not a multiple of shorter object length")
    } else  B <- rep(B,pdim %/% length(B))
  }
  if (!is.double(B[1]))  B <- as.double(B)
  return(array( .Fortran("subfullsparse",
                          nrow,ncol,dcheck(A@entries),A@colindices,
                          A@rowpointers,b=dcheck(B),PACKAGE = "spam"
                          )$b,c(nrow,ncol)))
}

".spam.subsparsefull" <- function(B,A){
  # A is sparse, B is full
  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  pdim <- prod(nrow,ncol)
  if (is.matrix(B)) {
    if(ncol != dim(B)[2] || nrow != dim(B)[1])
      stop("non-conformable matrices")
  } else {
    if(pdim %% length(B)!=0) {
      stop("longer object length
        is not a multiple of shorter object length")
    } else  B <- rep(B,pdim %/% length(B))
  }
  if (!is.double(B[1]))  B <- as.double(B)
  return(array( .Fortran("subsparsefull",
                          nrow,dcheck(A@entries),A@colindices,
                          A@rowpointers,b=dcheck(B),PACKAGE = "spam"
                          )$b,c(nrow,ncol)))
}

".spam.addsubsparsesparse" <-
function(A,B,s)
{
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  if(ncol != B@dimension[2] || nrow != B@dimension[1])
    stop("non-conformable matrices")
  ###IMPROVEME : there is a fortran routine getting the correct nr!

#      subroutine aplbdg (nrow,ncol,ja,ia,jb,ib,ndegr,nnz,iw) 
#      integer ja(*),jb(*),ia(nrow+1),ib(nrow+1),iw(ncol),ndegr(nrow) 
  
  nzmax <- .Fortran("aplbdg",
                nrow,                ncol,
                A@colindices,               A@rowpointers,
                B@colindices,               B@rowpointers,
                vector("integer",nrow),nnz=vector("integer",1),vector("integer",ncol),
                NAOK=!.Spam$safemode[3],DUP = FALSE,PACKAGE = "spam"
                )$nnz

  z <- .Fortran("aplsb1",
                nrow,
                ncol,
                dcheck(A@entries),
                A@colindices,
                A@rowpointers,
                as.double(s),
                dcheck(B@entries),
                B@colindices,
                B@rowpointers,
                entries     = vector("double",nzmax),
                colindices  = vector("integer",nzmax),
                rowpointers = vector("integer",nrow+1),
                as.integer(nzmax+1),
                ierr = vector("integer",1),
                NAOK=!.Spam$safemode[3],DUP = FALSE,PACKAGE = "spam"
                )
  if(z$ierr != 0) stop("insufficient space for sparse matrix addition")
  nz <- z$rowpointers[nrow+1]-1
  newz <- new("spam")
  slot(newz,"entries",check=FALSE) <- z$entries[1:nz]
  slot(newz,"colindices",check=FALSE) <- z$colindices[1:nz]
  slot(newz,"rowpointers",check=FALSE) <- z$rowpointers
  slot(newz,"dimension",check=FALSE) <- c(nrow,ncol)
  return(newz)
}


".spam.elemul" <-
function(e1,e2)
{
  if(is.vector(e1)) {
    if(length(e1) == 1){
      if(e1==0) return( spam(0,nrow(e2),ncol(e2)))
      else{  # just a scalar
        e2@entries <- e1*e2@entries
        return(e2)
      }
    }  else if(length(e1) == nrow(e2))
      return(diag.spam(e1) %*% e2)
    else # length(e1) == ncol(e2) is not required
      stop("e1 and e2 not conformable for efficient element-by-element multiplication")
  }
  else if(is.vector(e2)) {
    if(length(e2) == 1){
      if(e2==0)   return( spam(0,nrow(e1),ncol(e1)))
      else {
        e1@entries <- e2*e1@entries
        return(e1)
      }
    }
    else if(length(e2) == nrow(e1))
      return(diag.spam(e2) %*% e1)
    else
      stop("e1 and e2 not conformable for efficient element-by-element multiplication")
  }
  if(is.matrix(e1))
    e1 <- as.spam(e1)
  else if(is.matrix(e2))
    e2 <- as.spam(e2)
  if(!(is.spam(e1) && is.spam(e2)))
    stop("Arguments must be of class:  vector, matrix or spam")
  
  e1row <- e1@dimension[1]
  e1col <- e1@dimension[2]
  if(e1col != e2@dimension[2] | e1row != e2@dimension[1])
    stop("non-conformable matrices")
  nnzmax <- length(intersect(e1@colindices+e1col*(rep(1:e1row,diff(e1@rowpointers))-1),
                             e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1)))+1
  z <- .Fortran("aemub",
                e1row,
                e1col,
                dcheck(e1@entries),
                e1@colindices,
                e1@rowpointers,
                dcheck(e2@entries),
                e2@colindices,
                e2@rowpointers,
                entries     = vector("double",nnzmax),
                colindices  = vector("integer",nnzmax),
                rowpointers = vector("integer",e1row+1),
                integer(e1col),
                double(e1col),
                as.integer(nnzmax),
                ierr = vector("integer",1),
                DUP = FALSE,
                PACKAGE = "spam"
                )
  if(z$ierr != 0)      stop("insufficient space for element-wise sparse matrix multiplication")
  nnz <- z$rowpointers[e1row+1]-1
  if(identical(z$entries,0)){#trap zero matrix
    z$colindices <- as.integer(1)
    z$rowpointers <- as.integer(c(1,rep(2,nrow)))
  }
  return(new("spam",entries=z$entries[1:nnz],colindices=z$colindices[1:nnz],rowpointers=z$rowpointers,dimension=c(e1row,e1col)))
}

                                        
".spam.elediv" <-
function(e1,e2)
{ # Element-wise matrix division of two spams
if(is.numeric(e1) && length(e1) == 1)
  { e2@entries <- e1/e2@entries
    return(e2)
  } else if(is.numeric(e2) && length(e2) == 1) {
    e1@entries <- e1@entries/e2
    return(e1)
  }
  else if(is.spam(e1) || is.spam(e2) || is.matrix(e1) || is.matrix(e2)){
        if(is.matrix(e1)) e1 <- as.spam(e1)
        if(is.matrix(e2)) e2 <- as.spam(e2)
        nrow <- e1@dimension[1]
        ncol <- e1@dimension[2]
        if(ncol != e2@dimension[2] | nrow != e2@dimension[1])
          stop("matrices not conformable for element-by-element division")
	nzmax <- length(unique(c(e1@colindices+ncol*(rep(1:nrow,diff(e1@rowpointers))-1),
                                 e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1))))+1
        z <- .Fortran("aedib",
                      nrow,
                      ncol,
                      icheck(1),
                      dcheck(e1@entries),
                      e1@colindices,
                      e1@rowpointers,
                      dcheck(e2@entries),
                      e2@colindices,
                      e2@rowpointers,
                      entries     = vector("double",nzmax),
                      colindices  = vector("integer",nzmax),
                      rowpointers = vector("integer",nrow+1),
                      as.integer(nzmax),
                      integer(ncol),
                      double(ncol),
                      ierr = vector("integer",1),
                      DUP = FALSE,
                      PACKAGE = "spam"
                      )
        if(z$ierr != 0) stop("insufficient space for element-wise sparse matrix division")
        nz <- z$rowpointers[nrow+1]-1
        return(new("spam",entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,dimension=c(nrow,ncol)))
      }
  else stop("Arguments have to be class 'spam' or numeric")
}

                                        
".spam.expo" <-
  function(e1,e2)
{
# Performs element-wise exponentiation on sparse matrices
  if(is.numeric(e1) && length(e1) == 1) {
    e2@entries=e1^e2@entries
    return(e2)
  } else if(is.numeric(e2) && length(e2) == 1) {
    e1@entries=e1@entries^e2
    return(e1)
  } else if(is.spam(e1) || is.spam(e2) || is.matrix(e1) || is.matrix(e2)){
        if(is.matrix(e1)) e1 <- as.spam(e1)
        if(is.matrix(e2)) e2 <- as.spam(e2)
        nrow <- e1@dimension[1]
        ncol <- e1@dimension[2]
        if(ncol != e2@dimension[2] | nrow != e2@dimension[1])
          stop("matrices not conformable for element-wise exponentiation ")
	nzmax <- length(unique(c(e1@colindices+ncol*(rep(1:nrow,diff(e1@rowpointers))-1),
                                 e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1))))+1
        z <- .Fortran("aeexpb",
                icheck(nrow),
                icheck(ncol),
		as.integer(1),
                dcheck(e1@entries),
                icheck(e1@colindices),
                icheck(e1@rowpointers),
                dcheck(e2@entries),
                icheck(e2@colindices),
                icheck(e2@rowpointers),
                entries     = vector("double",nzmax),
                colindices  = vector("integer",nzmax),
                rowpointers = vector("integer",nrow+1),
                as.integer(nzmax),
		integer(ncol),
		double(ncol),
                ierr = vector("integer",1),
                PACKAGE = "spam"
                      )
        if(z$ierr != 0) stop("insufficient space for element-wise exponentiation")
        nz <- z$rowpointers[nrow+1]-1
        return(new("spam",entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,dimension=c(nrow,ncol)))
      }
  else stop("Arguments have to be class 'spam' or numeric")
}


# SUBSETTING
##########################################################################################

# notice the drop catch...
#   I don't know the best and official way, but it works as it is here...

setMethod("[", signature(x = "spam",
			 i = "missing", j = "missing", drop = "ANY"),
	  function (x, i, j, drop) { # cat("missmiss")
           x})

setMethod("[",signature(x="spam",i="vector",j="missing", drop = "ANY"),
	  function (x, i, j, drop) { # cat("vecmiss")
            subset.spam(x,rw=i,cl=1:x@dimension[2],drop)})

setMethod("[",signature(x="spam",i="vector",j="vector", drop = "ANY"),
	  function (x, i, j,drop) { # cat("vecvec")
            subset.spam(x,rw=i,cl=j,drop)} )

setMethod("[",signature(x="spam",i="missing",j="vector", drop = "ANY"),
	  function (x, i, j,drop) { # cat("missvec")
            subset.spam(x,rw=1:x@dimension[1],cl=j,drop)} )

setMethod("[",signature(x="spam",i="matrix",j="missing", drop = "ANY"),
	  function (x, i, j, drop) {subset.spam(x,rw=i,drop) })

setMethod("[",signature(x="spam",i="matrix",j="matrix", drop = "ANY"),
	  function (x, i, j, drop) {subset.spam(x,rw=cbind(c(i),c(j)),drop) })

setMethod("[",signature(x="spam",i="spam",j="missing", drop = "ANY"),
	  function (x, i, j, drop=.Spam$drop) 
{
  # drop is not implemented yet
  dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]
  if ( i@dimension[1]>nrow | i@dimension[2]>ncol)
    stop("subscript out of bounds",call.=FALSE)
  z <- .Fortran("amask",
                nrow=nrow,
                ncol=ncol,
                a=dcheck(x@entries),
                colindices=icheck(x@colindices),
                rowpointers=icheck(x@rowpointers),
                jmask=i@colindices,
                imask=c(i@rowpointers,rep(i@rowpointers[length(i@rowpointers)],nrow+1-length(i@rowpointers))),
                c=dcheck(x@entries),
                jc=icheck(x@colindices),
                ic=icheck(x@rowpointers),           
                iw=logical(ncol),
                nzmax=length(i@colindices) ,
                ierr=as.integer(0),
                PACKAGE="spam")
  nz <- z$ic[nrow+1]-1
  if (nz==0) return( numeric(0))
  if (drop) {
    ic <- unique( z$ic[1:(z$nrow+1)])
    dimx <- as.integer(c(length(ic)-1,max(z$jc[1:nz])))    
  } else {
    ic <-z$ic[1:(z$nrow+1)]
  }
  return(new("spam",entries=z$c[1:nz],colindices=z$jc[1:nz],rowpointers=ic,
               dimension=dimx))
}      )

setMethod("[", signature(x = "spam", i = "ANY", j = "ANY", drop = "ANY"),
	  function(x,i,j, drop)
          stop("Invalid or not-yet-implemented 'spam' subsetting"))

# the proper S3 subsetting causes problems... 
"[.spam" <- function (x, rw, cl,drop=.Spam$drop) {subset.spam(x,rw=rw,cl=cl,drop) }

"subset.spam" <-
function (x,rw,cl,drop=.Spam$drop,...)
{
  # we separate into cases where:
  # (A) rw matrix:
  #     1: logical: transformation to spam and extract structure
  #     2: two column matrix: extract (i,j) as given by the lines.
  #     3: all else extract   x[ c( rw)]
  # (B) rw and cl one element: ((i,j)
  # (C) rw and cl vectors:  (i1:i2,j1:j2)               [i1<=i2, j1<=j2]
  #                         (c(i1,...,ii),c(j1,...,jj)) [arbitrary block]

  if (missing(drop)) drop <- .Spam$drop
  dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]
  
  if (is.matrix(rw)) {
    if (is.logical(rw)) {
      return( x[as.spam.matrix(rw)] )
    }
    if (dim(rw)[2]==2) {
      ir <- rw[,1]
      jr <- rw[,2]
    } else  {
      ir <- c(rw-1) %% nrow + 1
      jr <- c(rw-1) %/% nrow + 1
    }
    if ( (min(ir)<1)|(max(ir)>x@dimension[1])|(min(jr)<1)|(max(jr)>x@dimension[2]))
      stop("subscript out of bounds",call.=FALSE)
    nir <- length(ir)
    return(.Fortran("getallelem",
                    nir,
                    as.integer(ir),
                    as.integer(jr),
                    dcheck(x@entries),icheck(x@colindices),icheck(x@rowpointers),
                    integer(nir),
                    allelem=vector("double",nir),
                    DUP=FALSE,
                      PACKAGE="spam")$allelem)

  }
  # negative values:
  if ( max(rw)<0 )       rw <- seq_len( nrow)[rw] 
  if ( max(cl)<0 )       cl <- seq_len( ncol)[cl] 
  
  # logical
  if (is.logical(rw))    rw <- seq_len( nrow)[rw] 
  if (is.logical(cl))    cl <- seq_len( ncol)[cl] 
  
  if (length(cl)==0) stop("You should subset at least one element for the columns",call.=FALSE)
  if (length(rw)==0) stop("You should subset at least one element for the rows",call.=FALSE)

  if ( (min(rw)<1)|(max(rw)>x@dimension[1])|(min(cl)<1)|(max(cl)>x@dimension[2]))
    stop("subscript out of bounds",call.=FALSE)
  
  if (length(rw)==1 & length(cl)==1){
                                        # function to extract only one element
    return(.Fortran("getelem",
                    as.integer(rw),
                    as.integer(cl),
                    dcheck(x@entries),icheck(x@colindices),icheck(x@rowpointers),
                    iadd=as.integer(0),
                    elem=as.double(0),
                    PACKAGE="spam")$elem)
  }
  if (is.vector(rw) && is.vector(cl)) {
    nrw <- length(rw)   # length returns an integer, so is a product therof
    ncl <- length(cl)
    diffrw <- diff(rw)
    diffcl <- diff(cl)
    nz <- as.integer( min( (1+sum(diff(sort(rw))==0))*(1+sum(diff(sort(cl))==0))*
                          length(x@entries), prod(nrw,ncl)))  # very pessimistic
    if (all(diffrw==1) & all(diffcl==1)) {
      z <- .Fortran("submat",
                    nrow,
                    job=as.integer(1), # need values as well
                    i1=as.integer(rw[1]),
                    i2=as.integer(rw[nrw]),
                    j1=as.integer(cl[1]),
                    j2=as.integer(cl[ncl]),
                    dcheck(x@entries),x@colindices,x@rowpointers,
                    nr=as.integer(0),
                    nc=as.integer(0),
                    entries=vector("double",nz),
                    colindices=vector("integer",nz),rowpointers=vector("integer",nrw+1),
                    NAOK=!.Spam$safemode[3],DUP = FALSE,PACKAGE = "spam")
      nz <- z$rowpointers[z$nr+1]-1
    } else {
      z <- .Fortran("getblock",
                    dcheck(x@entries),x@colindices,x@rowpointers,
                    nr=nrw,as.integer(rw),
                    nc=ncl,as.integer(cl),
                    nz=nz, entries=vector("double",nz),
                    colindices=vector("integer",nz),rowpointers=vector("integer",nrw+1),
                    NAOK=!.Spam$safemode[3],DUP = FALSE,PACKAGE = "spam")
      nz <- z$nz
    }
    if (nz==0) {#trap zero matrix
      if (drop==TRUE && (z$nr==1 || z$nc==1)) return( vector("double",max(z$nr,z$nc)))
      else
        return(new("spam",rowpointers=c(int1,rep.int(int2,z$nr )),
                   dimension = c(z$nr,z$nc)))
    }  
    
    if (drop==TRUE && (z$nr==1 || z$nc==1))
      # this is essentially a c() call
      return(.Fortran("spamcsrdns",
                 nrow=z$nr,
                 entries=z$entries[1:nz],
                 colindices=z$colindices[1:nz],
                 rowpointers=z$rowpointers[1:(z$nr+1)],
                 res=vector("double",prod(z$nr,z$nc)),  
                 NAOK=!.Spam$safemode[3],DUP=FALSE,PACKAGE = "spam")$res)
    else {
      newx <- new("spam")
      slot(newx,"entries",check=FALSE) <- z$entries[1:nz]
      slot(newx,"colindices",check=FALSE) <- z$colindices[1:nz]
      slot(newx,"rowpointers",check=FALSE) <- z$rowpointers[1:(z$nr+1)]
      slot(newx,"dimension",check=FALSE) <- c(z$nr,z$nc)
      return(newx)
    }
  
  }
  stop("invalid or not-yet-implemented 'spam' subsetting")
}

# ASSIGNING:
##########################################################################################

# as S3subsetting causes problems, we eliminate this... 
#"[<-.spam" <- function (x, rw, cl,value) {#cat('qq');
#                                          assign.spam(x,rw,cl,value) }


setReplaceMethod("[", signature(x = "spam",
			 i = "missing", j = "missing", value = "ANY"),
	  function (x, i, j, value) {#cat("mm");
                                     assign.spam(x,1:x@dimension[1],1:x@dimension[2],value)})

setMethod("[<-",signature(x="spam",i="vector",j="missing", value = "ANY"),
	  function (x, i, j, value) {#cat(i);
            assign.spam(x,i,1:x@dimension[2],value)} )

setMethod("[<-",signature(x="spam",i="vector",j="vector", value = "ANY"),
	  function (x, i, j, value) {#cat("vv");
            assign.spam(x,i,j,value)} )

setMethod("[<-",signature(x="spam",i="missing",j="vector", value = "ANY"),
	  function (x, i, j, value) {#cat(j);
            assign.spam(x,1:x@dimension[1],j,value)} )

setMethod("[<-",signature(x="spam",i="matrix",j="missing", value = "ANY"),
	  function (x, i, j, value) {#cat("Mm");
            assign.spam(x,i,NULL,value) })

setMethod("[<-",signature(x="spam",i="matrix",j="matrix",value = "ANY"),
	  function (x, i, j, value) {#cat("MM");
            assign.spam(x,cbind(c(i),c(j)),NULL,value) })

setMethod("[<-",signature(x="spam",i="spam",j="missing", value = "ANY"),
	  function (x, i, j, value) 
{
# cat("spam");
  dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]
  if ( i@dimension[1]>nrow | i@dimension[2]>ncol)
    stop("subscript out of bounds",call.=FALSE)
  if ( ( (i@rowpointers[i@dimension[1]+1]-1) %%length(value))!= 0)
      stop("number of items to replace is not a multiple of replacement length")
  nzmax <- as.integer(min(prod(nrow,ncol), i@rowpointers[i@dimension[1]+1]+x@rowpointers[nrow+1]-2))
  if (length(value)!=  (i@rowpointers[i@dimension[1]+1]-1) )
    value <- rep(value, (i@rowpointers[i@dimension[1]+1]-1) %/%length(value))
#   cat(length(value))#@@#
  z <- .Fortran("subass",
                nrow,ncol,
                dcheck(x@entries),      x@colindices,    x@rowpointers,
                b=as.double(value),  bj=i@colindices, bi=i@rowpointers,
                c=vector("double",nzmax),jc=vector("integer",nzmax),ic=vector("integer",nrow+1),
                nzmax=nzmax,
                PACKAGE="spam")
    cnz <- z$ic[nrow+1]-1
  return(new("spam",entries=z$c[1:cnz],colindices=z$jc[1:cnz],
             rowpointers=z$ic[1:(nrow+1)],dimension=c(nrow,ncol)))
}      )

setMethod("[<-", signature(x = "spam", i = "ANY", j = "ANY", value = "ANY"),
	  function(x,i,j, value){#  cat(value,class(value))
          stop("invalid or not-yet-implemented 'spam' subassigning")})



"assign.spam" <-
function (x, rw, cl,value)
{
  # we separate into cases where:
  # (A) rw matrix:
  #     1: logical: transformation to spam and extract structure
  #     2: two column matrix: extract (i,j) as given by the lines.
  #     3: all else extract   x[ c( rw)]
  # (B) rw and cl one element: ((i,j)
  # (C) rw and cl vectors:  (i1:i2,j1:j2)               [i1<=i2, j1<=j2]
  #                         (c(i1,...,ii),c(j1,...,jj)) [arbitrary block]

#  print(rw)
#  print(cl)
#  print(value)
  if (!is.numeric(value)) stop("Assignment of numeric structures only")
  
  dimx <- x@dimension
  nrow <- dimx[1]
  ncol <- dimx[2]
  
  if (is.matrix(rw)) {
    if (is.logical(rw)) {
      return( x[as.spam(rw)] <- value)
    }
    if (dim(rw)[2]==2) {
      ir <- rw[,1]
      jr <- rw[,2]
    } else  {
      ir <- c(rw-1) %% nrow + 1
      jr <- c(rw-1) %/% nrow + 1
      rw <- cbind(ir,jr)
    }
    if ( (min(ir)<1)|(max(ir)>x@dimension[1])|(min(jr)<1)|(max(jr)>x@dimension[2]))
      stop("subscript out of bounds",call.=FALSE)
    if (any(duplicated(cbind(ir,jr))))
      stop("only unique index for subassigning",call.=FALSE)
    nir <- length(ir)
    if (length(value)!=nir)
      stop("number of items to replace is not a multiple of replacement length")
    
    ord <- order(ir,jr)
    rw <- rw[ord,,drop=F]
    bia <- .Fortran("constructia",
                    nrow,as.integer(nir),
                    rowpointers=vector("integer",nrow+1),
                    ir=as.integer(c(rw[,1],0)),
                    PACKAGE="spam")$rowpointers
    nzmax <- as.integer(min(prod(nrow,ncol), nir+x@rowpointers[nrow+1])+2)
    z <- .Fortran("subass",
                  nrow,ncol,
                  dcheck(x@entries), x@colindices, x@rowpointers,
                  b=as.vector(value[ord],"double"),
                  bj=as.vector(rw[,2],"integer"),  bi=bia,
                  entries=vector("double",nzmax),
                  colindices=vector("integer",nzmax),
                  rowpointers=vector("integer",nrow+1),
                  nzmax=nzmax,
                  DUP=FALSE,
                  PACKAGE="spam")
    cnz <- z$rowpointers[nrow+1]-1
    if (cnz<0) {
      cat('Negative cnz in subassigning, forced to one. Please report.')
      return( spam(0))
    }
    newx <- new("spam")
    slot(newx,"entries",check=FALSE) <- z$entries[1:cnz]
    slot(newx,"colindices",check=FALSE) <- z$colindices[1:cnz]
    slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
    slot(newx,"dimension",check=FALSE) <- c(nrow,ncol)
    return(newx)
    
  }
  # negative subsetting:
  if ( max(rw)<0 )    rw <- seq_len( nrow)[rw] 
  if ( max(cl)<0 )    cl <- seq_len( ncol)[cl] 
  
  # logical
  if (is.logical(rw))    rw <- seq_len( nrow)[rw] 
  if (is.logical(cl))    cl <- seq_len( ncol)[cl] 

  # sanity check
  if (length(rw)==0) stop("You should assign at least one element for the rows",call.=FALSE)
  if (length(cl)==0) stop("You should assign at least one element for the columns",call.=FALSE)


  if ( (min(rw)<1)|(max(rw)>x@dimension[1])|(min(cl)<1)|(max(cl)>x@dimension[2]))
    stop("subscript out of bounds",call.=FALSE)
  
  if (is.vector(rw) && is.vector(cl)) {
    if (any(duplicated(rw))||any(duplicated(cl)))
      stop("only unique index for subassigning",call.=FALSE)

    nrw <- length(rw)   # length returns an integer, so is a product therof
    ncl <- length(cl)
    bnz <- nrw*ncl

    if ( (bnz%%length(value))!= 0)
      stop("number of items to replace is not a multiple of replacement length")

    # we pack the value into a vector _row by row_
    value <- c(t(array(as.double(value),c(nrw,ncl))[order(rw),order(cl)]))
    
    bia <- vector("integer",nrow)  # bia has size of nrow + 1
    bia[rw] <- ncl        # in each row we have ncl new objects
    bia <- as.integer(c(1,cumsum(bia)+1))
                
    # we construct now a sparse matrix containing the "value" at positions rw and cl.
    # then we use the subassign function.
    nzmax <- as.integer(min(prod(nrow,ncol), bnz+x@rowpointers[nrow+1])+2)
    # new("spam",entries=value,colindices=rep(sort(as.integer(cl)),nrw),rowpointers=bia,c(nrow,ncol))
    z <- .Fortran("subass",
                  nrow,ncol,
                  dcheck(x@entries), x@colindices ,x@rowpointers,
                  b=value,
                  bj=rep(sort(as.integer(cl)),nrw),
                  bi=bia,
                  entries=vector("double",nzmax),colindices=vector("integer",nzmax),
                  rowpointers=vector("integer",nrow+1),
                  nzmax=nzmax,
                  PACKAGE="spam")
    cnz <- z$rowpointers[nrow+1]-1
    newx <- new("spam")
    slot(newx,"entries",check=FALSE) <- z$entries[1:cnz]
    slot(newx,"colindices",check=FALSE) <- z$colindices[1:cnz]
    slot(newx,"rowpointers",check=FALSE) <- z$rowpointers
    slot(newx,"dimension",check=FALSE) <- c(nrow,ncol)
    return(newx)
  }
  stop("invalid or not-yet-implemented 'spam' subsetting")
}

".spam.matmul.mat" <-
function(x,y)
{
    nrow <- x@dimension[1]
    ncol <- x@dimension[2]
    yrow <- dim(y)[1]
    ycol <- dim(y)[2]
    if(yrow != ncol)stop("not conformable for multiplication")
    z <- .Fortran("amuxmat",
                  nrow,
		  yrow,
		  ycol,
                  as.double(y),
                  y=vector("double",nrow*ycol),
                  dcheck(x@entries),
                  x@colindices,
                  x@rowpointers,
                  NAOK=!.Spam$safemode[3],
                  DUP=FALSE,
                  PACKAGE = "spam")$y
    dim(z) <- c(nrow,ycol)
    return(z)
  }



".spam.matmul" <-
function(x,y)
{
  if (is.vector(x)) {
    y <- t(y)
    nrow <- y@dimension[1]
    ncol <- y@dimension[2]
    if(length(x) != ncol)  stop("not conformable for multiplication")
    z <- .Fortran("amux",
                  nrow,
                  as.double(x),
                  y=vector("double",nrow),
                  dcheck(y@entries),
                  y@colindices,
                  y@rowpointers,
                NAOK=!.Spam$safemode[3],
                DUP=FALSE,
                  PACKAGE = "spam")$y
    dim(z) <- c(1,nrow)
    return(z)
  } 
  if (is.vector(y)) {
    nrow <- x@dimension[1]
    ncol <- x@dimension[2]
    if(length(y) != ncol)stop("not conformable for multiplication")
    z <- .Fortran("amux",
                  nrow,
                  as.double(y),
                  y=vector("double",nrow),
                  dcheck(x@entries),
                  x@colindices,
                  x@rowpointers,
                NAOK=!.Spam$safemode[3],
                DUP=FALSE,
                  PACKAGE = "spam")$y
    dim(z) <- c(nrow,1)
    return(z)
  }
  if (is.matrix(y)) y <- as.spam(y)
  if (is.matrix(x)) x <- as.spam(x)


  #matrix multiply two sparse spam matrices

  xn <- x@dimension[1]
  xm <- x@dimension[2]
  yl <- y@dimension[2]
  if(xm != y@dimension[1])
    stop("matrices not conformable for multiplication")

  z <- .Fortran("amubdg",
                xn,xm,yl,
                x@colindices,x@rowpointers,
                y@colindices,y@rowpointers,
                integer(xn),
                nz = vector("integer",1),
                integer(yl),
                NAOK=!.Spam$safemode[3],
                DUP=FALSE,
                PACKAGE = "spam")
  nzmax <- z$nz
  z <- .Fortran("amub",
                xn,yl,
                as.integer(1),
                dcheck(x@entries), x@colindices, x@rowpointers,
                dcheck(y@entries), y@colindices, y@rowpointers,
                entries = vector("double",nzmax), colindices = vector("integer",nzmax),
                rowpointers = vector("integer",xn+1),
                as.integer(nzmax),
                integer(yl),
                ierr = vector("integer",1),
                NAOK=!.Spam$safemode[3],
                DUP=FALSE,
                PACKAGE = "spam")
  nz <- z$rowpointers[xn+1]-1
  if(z$ierr != 0) stop("insufficient space for sparse matrix multiplication")
  
      
  if(nz==0){#trap zero matrix
    z$entries <- 0
    z$colindices <- as.integer(1)
    z$rowpointers <- as.integer(c(1,rep(2,xn)))
  }  else  z <- .Fortran("sortrows",
                         xn,entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,
                         NAOK=!.Spam$safemode[3],
                         PACKAGE = "spam")
  newz <- new("spam")
  slot(newz,"entries",check=FALSE) <- z$entries
  slot(newz,"colindices",check=FALSE) <- z$colindices[1:nz]
  slot(newz,"rowpointers",check=FALSE) <- z$rowpointers
  slot(newz,"dimension",check=FALSE) <- icheck(c(xn,yl))
  return(newz)
}


if(getRversion() < "2.6") {


  
# The following is taken form https://svn.r-project.org/R-packages/trunk/Matrix/R/AllGeneric.R
###---- Group Generics ----
## The following are **WORKAROUND** s currently needed for all non-Primitives:

## [The following is more future-proof than direct  setGeneric(.) calls:
## FIX (in R!) : "trunc" should really be in Math, but we try both for the time

# the following covers all from 'getGroupMembers("Math")'
for(fname in intersect(getGroupMembers("Math"),
		       c("log", "log2", "log10", "logb", "log1p", "expm1",
			 "gamma", "lgamma", "digamma", "trigamma",
			 "cummax", "cummin", "trunc")))
    if(!is.primitive(get(fname))) setGeneric(fname, group="Math")

## "Math2"
for(fname in intersect(getGroupMembers("Math2"),
		       c("round", "signif", "trunc")))
    if (!is.primitive(get(fname))) setGeneric(fname, group="Math2")

rm(fname)

# the following covers all from 'getGroupMembers("Summary")
# for R2.3.x this should work.
# Martin proposes a workaround:
#     http://tolstoy.newcastle.edu.au/R/help/05/12/18192.html
#

  .max_def <- function(x, ..., na.rm = FALSE) base::max(x, ..., na.rm = na.rm)
  .min_def <- function(x, ..., na.rm = FALSE) base::min(x, ..., na.rm = na.rm)
  .range_def <- function(x, ..., na.rm = FALSE) base::range(x, ..., na.rm = na.rm)
  .prod_def <- function(x, ..., na.rm = FALSE) base::prod(x, ..., na.rm = na.rm)
  .sum_def <- function(x, ..., na.rm = FALSE) base::sum(x, ..., na.rm = na.rm)
  .any_def <- function(x, ..., na.rm = FALSE) base::any(x, ..., na.rm = na.rm)
  .all_def <- function(x, ..., na.rm = FALSE) base::all(x, ..., na.rm = na.rm)

  setGeneric("max", function(x, ..., na.rm = FALSE) standardGeneric("max"),
             useAsDefault = .max_def, group = "Summary")
  setGeneric("min", function(x, ..., na.rm = FALSE) standardGeneric("min"),
             useAsDefault = .min_def, group="Summary")
  setGeneric("range", function(x, ..., na.rm = FALSE) standardGeneric("range"),
             useAsDefault = .range_def, group="Summary")
  setGeneric("prod", function(x, ..., na.rm = FALSE) standardGeneric("prod"),
             useAsDefault = .prod_def, group="Summary")
  setGeneric("sum", function(x, ..., na.rm = FALSE) standardGeneric("sum"),
             useAsDefault = .sum_def, group="Summary")
  setGeneric("any", function(x, ..., na.rm = FALSE) standardGeneric("any"),
             useAsDefault = .any_def, group="Summary")
  setGeneric("all", function(x, ..., na.rm = FALSE) standardGeneric("all"),
             useAsDefault = .all_def, group="Summary")
}


setMethod("Math","spam", function(x){ x@entries <- callGeneric(x@entries);x })
setMethod("Math2",signature(x = "spam", digits = "ANY"),
          function(x, digits){ x@entries <- callGeneric(x@entries, digits = digits);x })

setMethod("Summary","spam", function(x,...,na.rm=FALSE){ callGeneric(x@entries,...,na.rm=FALSE) })


setMethod("%*%",signature(x="spam",y="spam"),    .spam.matmul)
setMethod("%*%",signature(x="spam",y="matrix"),  .spam.matmul.mat)
setMethod("%*%",signature(x="spam",y="numeric"), .spam.matmul)
setMethod("%*%",signature(x="matrix",y="spam"),  .spam.matmul)
setMethod("%*%",signature(x="numeric",y="spam"), .spam.matmul)

setMethod("+",signature(e1="spam",e2="spam"),
          function(e1,e2){ .spam.addsubsparsesparse(e1,e2,1)})
setMethod("+",signature(e1="spam",   e2="ANY"),
          function(e1,e2){ .spam.addsparsefull(e1,e2)})
setMethod("+",signature(e1="ANY",   e2="spam"),
          function(e1,e2){ .spam.addsparsefull(e2,e1)})

setMethod("-",signature(e1="spam",   e2="ANY"),
          function(e1,e2){ .spam.subfullsparse(e1,e2)})
setMethod("-",signature(e1="ANY",   e2="spam"),
          function(e1,e2){ .spam.subsparsefull(e1,e2)})
setMethod("-",signature(e1="spam",e2="spam"),
          function(e1,e2){ .spam.addsubsparsesparse(e1,e2,-1)})


setMethod("*",signature(e1="spam",e2="spam"), .spam.elemul)
setMethod("*",signature(e1="spam", e2="ANY"), .spam.elemul)
setMethod("*",signature(e1="ANY", e2="spam"), .spam.elemul)

setMethod("/",signature(e1="spam",e2="spam"), .spam.elediv)
setMethod("/",signature(e1="spam", e2="ANY"), .spam.elediv)
setMethod("/",signature(e1="ANY", e2="spam"), .spam.elediv)

setMethod("&",signature(e1="spam",e2="spam"), 
          function(e1,e2){ z <- .spam.elemul(e1,e2);z@entries <- rep(1,length(z@colindices));z})
setMethod("&",signature(e1="spam",e2="ANY"), 
          function(e1,e2){ z <- .spam.elemul(e1,e2);z@entries <- rep(1,length(z@colindices));z})
setMethod("&",signature(e1="ANY",e2="spam"), 
          function(e1,e2){ z <- .spam.elemul(e1,e2);z@entries <- rep(1,length(z@colindices));z})

setMethod("|",signature(e1="spam",e2="spam"), 
          function(e1,e2){ z <- .spam.addsubsparsesparse(e1,e2,1);z@entries <- rep(1,length(z@colindices));z})
setMethod("|",signature(e1="spam",e2="ANY"), 
          function(e1,e2){ .spam.addsparsefull(e1,e2)!=0})
setMethod("|",signature(e1="ANY",e2="spam"), 
         function(e1,e2){ .spam.addsparsefull(e2,e1)!=0})


#setMethod("Arith",signature(e1="spam",e2="numeric"),
#          function(e1,e2) {e1@entries <- callGeneric(e1@entries,e2); e1} )
#setMethod("Arith",signature(e1="numeric",e2="spam"),
#          function(e1,e2) {e2@entries <- callGeneric(e2@entries,e1); e2} )
#setMethod("Arith",signature(e1="spam",e2="missing"),
#          function(e1,e2) {e1@entries <- callGeneric(e1@entries); e1} )

# The binary results are coerced to doubles!
# setMethod("Compare",signature(e1="spam",e2="numeric"),
#           function(e1,e2) {e1@entries <- as.double(callGeneric(e1@entries,e2)); e1} )
# setMethod("Compare",signature(e1="numeric",e2="spam"),
#           function(e1,e2) {e2@entries <- as.double(callGeneric(e2@entries,e1)); e2} )

#####################################################################################

upper.tri.spam <- function(x,diag=FALSE)
  {
    dimx <- x@dimension
    nrow <- dimx[1]
    z <- .Fortran("getu",
                  nrow,
                  dcheck(x@entries),x@colindices,x@rowpointers,
                  entries=dcheck(x@entries),colindices=x@colindices,rowpointers=x@rowpointers,
                PACKAGE="spam")
    nz <- z$rowpointers[dimx[1]+1]-1
    if (!diag) {
      z <- .Fortran("getdia",
                      n=nrow,
                      m=nrow,
                      job=as.integer(1),
                      entries=z$entries[1:nz],
                      colindices=z$colindices[1:nz],
                      rowpointers=z$rowpointers,
                      len=nrow,
                      diag=vector("double",nrow),
                      idiag=vector("integer",nrow),
                      ioff=as.integer(0),
                PACKAGE = "spam"
                )
      nz <- z$rowpointers[nrow+1]-1
    }
    if(.Spam$trivalues)
      return(new("spam",entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,dimension=dimx))
    else
      return(new("spam",entries=rep(1,nz),colindices=z$colindices[1:nz],rowpointers=z$rowpointers,dimension=dimx))
  }

lower.tri.spam <- function(x,diag=FALSE)
{
  dimx <- x@dimension
  nrow <- dimx[1]
  z <- .Fortran("getl",
                nrow,
                dcheck(x@entries),x@colindices,x@rowpointers,
                entries=dcheck(x@entries),colindices=x@colindices,rowpointers=x@rowpointers,
                PACKAGE="spam")
  nz <- z$rowpointers[nrow+1]-1
  
  if (!diag) {
    z <- .Fortran("getdia",
                  n=nrow,
                  m=nrow,
                  job=as.integer(1),
                  entries=z$entries[1:nz],
                  colindices=z$colindices[1:nz],
                  rowpointers=z$rowpointers,
                  len=nrow,
                  diag=vector("double",nrow),
                  idiag=vector("integer",nrow),
                  ioff=as.integer(0),
                  PACKAGE = "spam"
                  )
    nz <- z$rowpointers[nrow+1]-1
  }
  if(.Spam$trivalues)
    return(new("spam",entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,dimension=dimx))
  else
    return(new("spam",entries=rep(1,nz),colindices=z$colindices[1:nz],rowpointers=z$rowpointers,dimension=dimx))
}



setGeneric("upper.tri")
setMethod("upper.tri","spam",upper.tri.spam)
setGeneric("lower.tri")
setMethod("lower.tri","spam",lower.tri.spam)

########################################################################

norm <- function(x, type = "sup", ...){
  typ <- charmatch(tolower(type), c("sup",'l1',"frobenius","hs"))
  if (is.na(typ))          stop("undefined norm '",type,"'.",call.=FALSE)

  switch(typ,
         max(abs(x)),
         sum(abs(x)),
         sqrt(sum(x^2)),sqrt(sum(x^2))
         )
}
norm.spam <- function(x, type = "sup", ...){
  typ <- charmatch(tolower(type), c("sup",'l1',"frobenius","hs"))
  if (is.na(typ))          stop("undefined norm '",type,"'.",call.=FALSE)

  switch(typ,
         max(abs(x@entries)),
         sum(abs(x@entries)),
         sqrt(sum(x@entries^2)),
         sqrt(sum(x@entries^2))
         )
}

setGeneric("norm",function(x, type = "sup",...)standardGeneric("norm"))
setMethod("norm",signature(x="spam",type="character"), norm.spam)
setMethod("norm",signature(x="spam",type="missing"),
          function(x,type) norm.spam(x,"sup"))




# fields uses the construct of vector representation for a diagonal matrix.

# Create a special matrix multiply for diagonal matrices.
# Diagonal matrix assumed to be just a vector.
# NOTE: this is not a symmetric operation:
#  when a left vector is given it is a diagonal matrix
#  when a right vector is given it is a vector.
#
.spam.diagmulmat <- function(x,y){
  nrow <- y@dimension[1]
  if(length(x) != nrow)  stop("not conformable for multiplication")
  z <- .Fortran("diagmua",
                nrow,
                entries=dcheck(y@entries),
                y@rowpointers,
                as.vector(x,"double"),
                PACKAGE = "spam")$entries
  y@entries <- z
  return(y)
}

.spam.diagaddmat <- function(x,y){
#      subroutine diagaddmat (nrow, a, ja, ia, diag, b, jb, ib, iw)
  nrow <- y@dimension[1]
  minrc <- min( y@dimension)
  if(length(x) != minrc)  stop("not conformable for addition")
  z <- .Fortran("diagaddmat",
                nrow = nrow,
                n = minrc,
                a = c(y@entries,double(minrc)),
                ja = c(y@colindices,integer(minrc)),
                ia = y@rowpointers,
                diag = as.double(x),
                iw = vector("integer",nrow),    # just to be sure
                info = vector("integer",nrow+1),
                PACKAGE = "spam")
  nz <- z$ia[nrow+1]-1
  return(new("spam",entries=z$a[1:nz],colindices=z$ja[1:nz],
             rowpointers=z$ia,dimension=y@dimension))
}

setGeneric("%d*%",function(x,y,...)standardGeneric("%d*%"))

setMethod("%d*%",signature(x="matrix",y="ANY"),       function(x,y){x%*%y} )
setMethod("%d*%",signature(x="numeric",y="matrix"),   function(x,y){x*y} )
setMethod("%d*%",signature(x="numeric",y="numeric"),  function(x,y){cbind(x*y)} )

setMethod("%d*%",signature(x="spam",y="spam"),    .spam.matmul )
setMethod("%d*%",signature(x="spam",y="ANY"),     .spam.matmul )
setMethod("%d*%",signature(x="numeric",y="spam"), .spam.diagmulmat )


setGeneric("%d+%",function(x,y,...)standardGeneric("%d+%"))
setMethod("%d+%",signature(x="matrix",y="ANY"),      function(x,y){ x+y } )
setMethod("%d+%",signature(x="numeric",y="matrix"),  function(x,y){ diag(x)+y} )
setMethod("%d+%",signature(x="numeric",y="numeric"), function(x,y){ diag(x)+y} )

setMethod("%d+%",signature(x="spam",y="spam"),     function(x,y){ .spam.addsubsparsesparse(e1,e2,1)})
setMethod("%d+%",signature(x="spam",y="ANY"),      function(x,y){ .spam.addsparsefull(e1,e2)})
setMethod("%d+%",signature(x="numeric",y="spam"),  .spam.diagaddmat )


#####################################################################
#
# a bit of matrix handling

all.equal.spam <- function (target, current, tolerance = .Machine$double.eps^0.5,
    scale = NULL, check.attributes = FALSE,...)
{
    if (check.attributes)
        warning("attributes are not supported for 'spam' objects. Ignoring 'check.attributes' argument")
    if (!is.spam(target)) stop("'target' should be of class 'spam'")    
    if (!is.spam(current)) {
        return(paste("target is spam, current is ", data.class(current), sep = ""))
    }
    msg <- NULL
    lt <- length(target)
    lc <- length(current)
    if (lt != lc) {
      return(paste("Lengths (", lt, ", ", lc, ") differ", sep = ""))
    }
    dt <- target@dimension
    dc <- current@dimension
    if ( !all( dt == dc ))
      return(paste("Dimensions ([",dt[1],",",dt[2],"], [",
                    dc[1],",",dc[2], "]) differ", sep = ""))
    tmp <- sum(target@colindices != current@colindices) 
    if ( tmp>0)
      msg <- c(msg,paste("Column-sparsity structure differ (at least",
                    tmp,"instance(s))"))
    
    tmp <- sum(target@rowpointers != current@rowpointers) 
    if ( tmp>0)
      msg <- c(msg,paste("Row-sparsity structure differ (at least",
                    tmp,"instance(s))"))

    xy <- mean(abs(target@entries - current@entries))
    what <- if (is.null(scale)) {
        xn <- mean(abs(target@entries))
        if (is.finite(xn) && xn > tolerance) {
            xy <- xy/xn
            "relative"
        }
        else "absolute"
    }
    else {
        xy <- xy/scale
        "scaled"
    }
    if (is.na(xy) || xy > tolerance)
        msg <- c(msg,paste("Mean", what, "difference:",
            format(xy)))
    if (is.null(msg))
        TRUE
    else msg

}

isSymmetric.spam <- function(object, tol = 100 * .Machine$double.eps, ...)
{
  # very similar to is.Symmetric.matrix
  if (!is.spam(object)) return(FALSE)
  
  d <- object@dimension
  if (d[1] != d[2])     return(FALSE)
  test <-  all.equal.spam(object, t.spam(object), tolerance = tol, ...)
  isTRUE(test)

}

setMethod("all.equal",signature(target="spam",current="spam"), all.equal.spam )
setMethod("all.equal",signature(target="matrix",current="spam"),
          function (target, current, tolerance = .Machine$double.eps^0.5,
                    scale = NULL, check.attributes = FALSE,eps = .Spam$eps,...)
{
    if (check.attributes)
      warning("attributes are not supported for 'spam' objects. Ignoring 'check.attributes' argument")
    msg <- NULL
    
    dimx <- dim(target)
    nz <- length(target)
    z <- .Fortran("spamdnscsr", nrow = dimx[1], ncol = dimx[2],
        x = as.double(target), dimx[1], entries = vector("double",
            nz), colindices = vector("integer", nz), rowpointers = vector("integer",
            dimx[1] + 1), eps = as.double(eps), NAOK = !.Spam$safemode[3],
        DUP = FALSE, PACKAGE = "spam")

    lt <- z$rowpointers[dimx[1] + 1] - 1
    lc <- length(current)
        
    if (lt != lc) {
      return(paste("Lengths (", lt, ", ", lc, ") differ", sep = ""))
    }
    dt <- dim(target)
    dc <- current@dimension
    if ( !all( dt == dc ))
      return(paste("Dimensions ([",dt[1],",",dt[2],"], [",
                    dc[1],",",dc[2], "]) differ", sep = ""))
    tmp <- sum(z$colindices[1:lt] != current@colindices) 
    if ( tmp>0)
      msg <- c(msg,paste("Column-sparsity structure differ (at least",
                    tmp,"instance(s))"))
    
    tmp <- sum(z$rowpointers != current@rowpointers) 
    if ( tmp>0)
      msg <- c(msg,paste("Row-sparsity structure differ (at least",
                    tmp,"instance(s))"))

    xy <- mean(abs(z$entries[1:lt] - current@entries))
    what <- if (is.null(scale)) {
        xn <- mean(abs(z$entries))
        if (is.finite(xn) && xn > tolerance) {
            xy <- xy/xn
            "relative"
        }
        else "absolute"
    }
    else {
        xy <- xy/scale
        "scaled"
    }
    if (is.na(xy) || xy > tolerance)
        msg <- c(msg,paste("Mean", what, "difference:",
            format(xy)))
    if (is.null(msg))
        TRUE
    else msg

}
 )
setMethod("isSymmetric","spam", isSymmetric.spam)
