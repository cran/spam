".onLoad" <-
function (lib, pkg)
{
  library.dynam("spam",pkg, lib)
# this is a workaround from the NAMESPACE inclusion  
assign(".Spam",list(eps=.Machine$double.eps,   # smaller than this is considered as zero
              drop=FALSE,                # drop passed to subset functions
              printsize=100,             # the max size which we print as regular matrices
              imagesize=10000,           # the max size which we display as regular matrices
              trivalues=FALSE,           # with upper./lower/.tri return values (TRUE) or only structure?
              cex=1200,                  # scaling factor for scatter displays

              version=list(major=0,
                minor=.11,
                year=2007,
                month=09,
                day=04),

              safemode=TRUE,             # verify double and integer formats. 
              bcksl=TRUE                 # what type of back/forwardsolve?
              
              ), env = .GlobalEnv)
}

".onAttach" <-
function (lib, pkg)
{
  cat("Package 'spam' is loaded.  Version ",
      .Spam$version$major+.Spam$version$minor," (",
      .Spam$version$year,"-",
      sprintf("%02d",.Spam$version$month),"-",
      sprintf("%02d",.Spam$version$day),").",
      "\nType demo( spam) for some demos,",
      " help( Spam) for an overview of this library.\n",
      sep='',fill=TRUE)

}
todo <- function() help( "spam.todo")
spam.history <- function() help("spam.history")




dcheck <- function(x) if (.Spam$safemode) as.double(x) else x
icheck <- function(x) if (.Spam$safemode) as.integer(x) else x

validspamobject <- function(object) {
  if(!(length(object@dimension) == 2) ){
    print(object@dimension)
    return("invalid dimension attribute")
  }
  else{
    nrow <- object@dimension[1]
    ncol <- object@dimension[2]
  }
  if (any(!is.double(object@entries)))
    return("matrix entries need to be of double mode")
  if (sum(is.infinite(object@entries))>0)
    return("'NA/NaN/Inf' not allowed")
  if(!(length(object@entries) ==length(object@colindices)))
    return("entries and column indices don't have equal lengths")
  if(any(object@colindices < 1) || any(object@colindices > ncol)) 
    return("column indices exceeds dimension bounds")
  if(any(object@rowpointers < 1))
    return("some elements of row pointes are <= 0")
  if(any(diff(object@rowpointers)<0))
    return("row pointers are not monotone increasing")
  diffcolindices <- diff(object@colindices)     # positive values within each row
  if (all(diff(object@rowpointers)>1) && length(diffcolindices)>0)   # only if we have multiple values
    if (nrow==1) {
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
  TRUE
}

setClass("spam",representation(entries="numeric",
                               colindices="integer",
                               rowpointers="integer",
                               dimension="integer"),
         validity = validspamobject)

setMethod("initialize", "spam", 
          function(.Object,
                   entries    = 0,         # by default a 1x1 zero matrix.
                   colindices = as.integer( rep(1, length( entries ))),  # or a nx1 matrix
                   rowpointers= as.integer( 1:(length( entries )+1)),   # with n=length(ra)
                   dimension  = as.integer( c(length( rowpointers )-1,max( colindices ))))
          {
            # a few specific "degenerate" cases: 
            if (length(entries)==0) entries <- 0         # e.g., induced by rep(1,0)
            if (rowpointers[ length(rowpointers)] ==1) {  # e.g., zero matrix
              rowpointers[-1] <- as.integer(2)
              colindices <- as.integer(1)
            }
            .Object@entries     <- entries
            .Object@colindices  <- colindices
            .Object@rowpointers <- rowpointers
            .Object@dimension   <- dimension
            validObject(.Object)
            .Object
          })


print.spam <- function(x,...) {
  if (prod(x@dimension) < .Spam$printsize) {
    print(as.matrix(x),...)
  } else {
    if ( (length(x@entries)==1)  & (x@entries[1]==0)) {
      cat(paste('Zero matrix of dimension ',x@dimension[1],'x',
                x@dimension[2],'.\n',sep=''))
    }else {
      cat(paste('Matrix of dimension ',x@dimension[1],'x',
                x@dimension[2],' with nonzero elements:\n',sep=''))
      print(x@entries,...)
    }
  }
  cat("Class 'spam'\n")
  invisible(x)
}

summary.spam <- function(object,...) {
            nz <- length(object@entries)
            dens <- nz/prod(object@dimension)*100
            cat("Matrix object of class 'spam' of ")
            cat(paste('dimension ',object@dimension[1],'x',object@dimension[2],'.\n',sep=''))
            cat(paste('  Density of the matrix is ',signif(dens,3),'% (nz=',nz,').\n',sep=''))
            cat("Class 'spam'\n")
            invisible(object)
          }


setMethod("show","spam",  function(object) {
  if (prod(object@dimension) < .Spam$printsize) {
    print(as.matrix(object))
  } else {
    if ( (length(object@entries)==1)  & (object@entries[1]==0)) {
      cat(paste('Zero matrix of dimension ',object@dimension[1],'x',
                object@dimension[2],'.\n',sep=''))
    }else {
      cat(paste('Matrix of dimension ',object@dimension[1],'x',
                object@dimension[2],' with nonzero elements:\n',sep=''))
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
setMethod("dim<-",   "spam",function(x,value) stop("operation not allowed on 'spam' object") )


setMethod("c","spam", function(x,...,recursive=TRUE){
  if (length( list(...)) < 1)
    c(as.matrix(x))
  else
    c(as.matrix(x),c(...,recursive),recursive)
})

########################################################################
# diag and derivatives
"diag.spam" <-
function(x=1, nrow, ncol=n)
{
  if (is.spam(x)) return(diag(x, nrow, ncol))

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
  p <- as.integer(ncol)

  m <- min(n, p)
  tmp <- numeric(m)
  tmp[1:m] <- as.double(x)
 
  return(new("spam",entries=tmp,   # numeric return a double
             colindices=1:m,rowpointers=as.integer(c(1:m,rep(m+1,n+1-m))),
             dimension=c(n,p)))
  
}



"diag.spam<-" <-
function(x,value)
{
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
                iw = integer(nrow),    # just to be sure
                info = integer(nrow+1),
                PACKAGE = "spam")
  nz <- z$ia[nrow+1]-1
  return(new("spam",entries=z$a[1:nz],colindices=z$ja[1:nz],
             rowpointers=z$ia,dimension=x@dimension))
}

"diag.of.spam" <-
function(x, nrow, ncol)
{
  len <- min(x@dimension)
  return(.Fortran("getdiag",
                a = dcheck(x@entries),
                colindices =  x@colindices,
                rowpointers = x@rowpointers,
                len = len,
                diag = double(len),
                PACKAGE = "spam"
                )$diag)
}

setMethod("diag","spam",diag.of.spam)
setMethod("diag<-","spam",get("diag.spam<-"))

########################################################################

"t.spam" <- function(x){
  dimx <- x@dimension
  nz <- x@rowpointers[dimx[1]+1]-1
  z <- .Fortran("smmptransp",
                n=dimx[1],
                m=dimx[2],
                ia=x@rowpointers,
                ja=x@colindices,
                diaga=as.integer(0),
                a=dcheck(x@entries),
                ib=integer(dimx[2]+1),
                jb=integer(nz),
                b=double(nz),
                move=as.integer(1),
                PACKAGE = "spam" 
                )
  return(new("spam",entries = z$b[1:nz], colindices = z$jb[1:nz], rowpointers = z$ib, dimension = dimx[2:1] ))
}
setMethod("t","spam",t.spam)

########################################################################

"is.spam" <- function(x) is(x,"spam")

"as.spam" <- function(x,  eps = .Spam$eps) stop('coercion not defined form this class')

"spam" <- function(x, nrow = 1, ncol = 1, eps = .Spam$eps) stop("argument 'x' should be of mode 'numeric' (or 'spam')")

"as.spam.spam" <-
function(x, eps = .Spam$eps)
{
  if (eps<0) stop("'eps' cannot be negative",call.=FALSE)
  dimx <- x@dimension
  z <- .Fortran("cleanspam",
                nrow=dimx[1],
                entries=dcheck(x@entries),
                colindices=x@colindices,
                rowpointers=x@rowpointers,
                eps=as.double(eps),
                PACKAGE = "spam"
                )
  nz <- z$rowpointers[dimx[1]+1]-1
  if(nz==0){#trap zero matrix
    z$entries <- 0
    z$colindices <- as.integer(1)
    z$rowpointers <- as.integer(c(1,rep(2,dimx[1])))
  }  
  return(new("spam", entries= z$entries[1:nz], colindices = z$colindices[1:nz],
             rowpointers = z$rowpointers[1:(dimx[1]+1)],    dimension = dimx))
}


"as.spam.matrix" <-
function(x, eps = .Spam$eps)
{  
  dimx <- dim(x)
  if (sum(is.infinite(x))>0) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[is.infinite(x)] <- 0
  }
  nz <- sum(abs(x)>eps)
  if(nz==0) return(new("spam",entries=0,colindices=as.integer(1),rowpointers=as.integer(c(1,rep(2,dimx[1]))), dimension=dimx))
                                        # no nonzero values. We preserve the dimension of x
  z <- .Fortran("spamdnscsr",
                nrow=dimx[1],
                ncol=dimx[2],
                nzmax=as.integer(nz),
                x=as.double(x),
                dimx[1],
                entries=double(nz),
                colindices=integer(nz),
                rowpointers=integer(dimx[1]+1),
                eps=as.double(eps),
                PACKAGE = "spam"
                )
  return(new("spam",entries = z$entries[1:nz], colindices = z$colindices[1:nz], rowpointers = z$rowpointers, dimension = dimx))
}

"as.spam.numeric" <-
function(x, eps = .Spam$eps)
{
  if (sum(is.infinite(x))>0) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[is.infinite(x)] <- 0
  }
  nz <- sum(abs(x)>eps)
  lx <- length(x)
  if(lx==1) # scalar
    return(new("spam",entries=as.double(ifelse(nz==0,0,x)), colindices = as.integer(1),
             rowpointers = as.integer(c(1,2)), dimension = as.integer(c(1,1))) ) 

  
  # standard vector.
  #   "initialize" handles only zero entries 
  return(new("spam",entries = as.double(x[abs(x)>eps]), colindices = rep(as.integer(1), nz),
             rowpointers = as.integer(cumsum(c(1,abs(x)>eps))), dimension = as.integer(c(lx,1))))
}



"spam.numeric" <-
function(x, nrow = 1, ncol = 1, eps = .Spam$eps)
{
  if (sum(is.infinite(x))>0) {
    warning("'NA/NaN/Inf' coerced to zero")
    x[is.infinite(x)] <- 0
  }
                                        # no modification thereof
  if (missing(nrow))
    nrow <- ceiling(length(x)/ncol)
  else if (missing(ncol))
    ncol <- ceiling(length(x)/nrow)
  if (length(x) == nrow * ncol)
    dim(x) <- c(nrow, ncol)
  else{
    if(length(x)==1 && abs(x)<eps) {
      dimx <- c(nrow,ncol)
      return(new("spam",entries=0,colindices=as.integer(1),
                 rowpointers=as.integer(c(1,rep(2,dimx[1]))), 
                 dimension=as.integer(dimx)))
    }
    else if((nrow*ncol)%%length(x)!=0)
      warning("ncol*nrow indivisable by length(x)")
      
    x <- rep(x, ceiling(nrow*ncol/length(x)))
    dim(x) <- c(nrow, ncol)
  }
  return( as.spam(x, eps=eps))
}

"spam.spam" <-
function(x, nrow = 1, ncol = 1, eps = .Spam$eps)
{
  x <- as.matrix(x)
  if (missing(nrow))
    nrow <- ceiling(length(x)/ncol)
  else if (missing(ncol))
    ncol <- ceiling(length(x)/nrow)
  if (length(x) == nrow * ncol)
    dim(x) <- c(nrow, ncol)
  else{
    if(length(x)==1 && abs(x)<eps) {
      dimx <- c(nrow,ncol)
      return(new("spam",entries=0,colindices=as.integer(1),
                 rowpointers=as.integer(c(1,rep(2,dimx[1]))), 
                 dimension=as.integer(dimx)))
    }
    else if((nrow*ncol)%%length(x)!=0)
      warning("ncol*nrow indivisable by length(x)")
      
    x <- rep(x, ceiling(nrow*ncol/length(x)))
    dim(x) <- c(nrow, ncol)
  }
  return( as.spam(x, eps=eps))
}



setGeneric("as.spam")
setMethod("as.spam","spam",   as.spam.spam)
setMethod("as.spam","matrix", as.spam.matrix)
setMethod("as.spam","numeric",as.spam.numeric)

setGeneric("spam")
setMethod("spam","numeric",spam.numeric)
setMethod("spam","spam",spam.spam)


########################################################################

if(paste(R.version$major, R.version$minor, sep=".") >= "2.5")
  {
    
"as.matrix.spam" <-
function(x,...){
  dimx <- x@dimension
  return(array(.Fortran("spamcsrdns",
                        nrow=dimx[1],
                        entries=dcheck(x@entries),
                        colindices=x@colindices,
                        rowpointers=x@rowpointers,
                        res=numeric(prod(dimx)),  # numeric is double! 
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
                        res=numeric(prod(dimx)),  # numeric is double! 
                        PACKAGE = "spam"
                        )$res,
               dimx)      # we preserve dimensions
         )
}
}


setMethod("as.matrix","spam",as.matrix.spam)

########################################################################

"rbind.spam" <-
function(...,deparse.level=0)
{
  if (deparse.level!=0) warning("Only 'deparse.level=0' implemented, coerced to zero,")
  addnargs <- ifelse(missing(deparse.level),0,1)

  nargs <- nargs()-addnargs
  if (nargs == 0)     return( NULL)
  args <- list(...)
  if (!is.null( names( args)))  {
    warning("Names of arguments are ignored")
    names( args) <- NULL
  }

  if (nargs == 1)     return( args[[1]])
  if (nargs == 2) {
    # this is the quick way
    if(!(is.spam(args[[1]])&is.spam(args[[2]])))
         stop("Not all argument are of class 'spam', in rbind.spam()",call.=FALSE)
    if(ncol(args[[1]])!=ncol(args[[2]]))
         stop("Arguments have differing numbers of columns, in rbind.spam()",call.=FALSE)
    
    nrow1 <- args[[1]]@dimension[1] 
    return(new("spam", entries = c(args[[1]]@entries, args[[2]]@entries),
               colindices =  c(args[[1]]@colindices,  args[[2]]@colindices),
               rowpointers = c(args[[1]]@rowpointers,
                 args[[2]]@rowpointers[-1]+args[[1]]@rowpointers[nrow1+1]-as.integer(1)),
               dimension =   c(nrow1+args[[2]]@dimension[1],args[[1]]@dimension[2])))
  } else {
    # "recursive" approach only, e.g. no checking
    tmp <- rbind.spam( args[[1]],args[[2]])
    for ( i in 3:nargs)
      tmp <- rbind.spam( tmp,args[[i]])
    return( tmp)
  }
}
  
  
"cbind.spam" <-
function(...,deparse.level=0)
{
  if (deparse.level!=0) warning("Only 'deparse.level=0' implemented, coerced to zero,")
  addnargs <- ifelse(missing(deparse.level),0,1)

  nargs <- nargs()-addnargs
  if (nargs == 0)     return( NULL)
  args <- list(...)
  if (!is.null( names( args)))  {
    warning("Names of arguments are ignored")
    names( args) <- NULL
  }

  if (nargs == 1)     return( args[[1]])
  if (nargs == 2) {
    # this is the quick way
    if(!(is.spam(args[[1]])&is.spam(args[[2]])))
         stop("Not all argument are of class 'spam', in cbind.spam()",call.=FALSE)
    if(nrow(args[[1]])!=nrow(args[[2]]))
         stop("Arguments have differing numbers of rows, in cbind.spam()",call.=FALSE)
    
    ncol1 <- args[[1]]@dimension[2] 
    nrow <-  args[[1]]@dimension[1] 

    rowpointers <- args[[1]]@rowpointers + (args[[2]]@rowpointers-as.integer(1))
    entries <- colindices <- NULL
    
    for (i in 1:nrow) {
      if (args[[1]]@rowpointers[i]<args[[1]]@rowpointers[i+1])
        stend1 <- args[[1]]@rowpointers[i]:(args[[1]]@rowpointers[i+1]-1)
      else stend1 <- NULL
      if (args[[2]]@rowpointers[i]<args[[2]]@rowpointers[i+1])
        stend2 <- args[[2]]@rowpointers[i]:(args[[2]]@rowpointers[i+1]-1)
      else stend2 <- NULL
      entries <- c( entries, args[[1]]@entries[stend1], args[[2]]@entries[stend2])
      colindices <-  c(  colindices, args[[1]]@colindices[stend1], args[[2]]@colindices[stend2]+ncol1)
    }
    return(new("spam", entries = entries,
               colindices =  colindices,
               rowpointers = rowpointers,
               dimension =   c(nrow, ncol1+args[[2]]@dimension[2])))
  } else {
    # "recursive" approach only, e.g. no checking
    tmp <- cbind.spam( args[[1]],args[[2]])
    for ( i in 3:nargs)
      tmp <- cbind.spam( tmp,args[[i]])
    return( tmp)
  }
}
  
  


setMethod("rbind","spam",rbind.spam)
setMethod("cbind","spam",cbind.spam)


########################################################################


".spam.compl" <- function(x){
# this function returns the structure of the zeros of the spam object x. 
  nrow <- x@dimension[1]
  ncol <- x@dimension[2]
  nnz <- x@rowpointers[nrow+1]-1
  nz <- nrow*ncol - nnz
  
  # we work through special cases      
  if(nz == 0)	            return(as.spam(0,nrow,ncol))
  if(nnz == 1 && x@entries == 0) return(as.spam(1,nrow,ncol))
  # normal case, proceed to efficient function
  z <- .Fortran("nzero",
                  dcheck(x@entries),
                  x@colindices,
                  x@rowpointers,
                  as.integer(nrow),
                  ncol,
                  as.integer(nnz),
                  as.integer(nz),
                  entries = double(nz),
                  colindices = integer(nz),
                  rowpointers = integer(nrow+1),
                  logical(ncol),
                PACKAGE="spam"
                )
    return(new("spam",entries=z$entries,colindices=z$colindices,rowpointers=z$rowpointers,dimension=x@dimension))
}

".spam.addsparsefull" <- function(A,B){
  # A is sparse, B is full
  if (missing(B)) return(A)
  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  pdim <- nrow*ncol
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
  return(matrix( .Fortran("addsparsefull",
                          nrow,dcheck(A@entries),A@colindices,
                          A@rowpointers,b=dcheck(B),PACKAGE = "spam"
                          )$b,nrow,ncol))
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
  pdim <- nrow*ncol
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
  return(matrix( .Fortran("subfullsparse",
                          nrow,ncol,dcheck(A@entries),A@colindices,
                          A@rowpointers,b=dcheck(B),PACKAGE = "spam"
                          )$b,nrow,ncol))
}
".spam.subsparsefull" <- function(B,A){
  # A is sparse, B is full
  if (!is.numeric(B)) stop("numeric argument expected")
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  pdim <- nrow*ncol
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
  return(matrix( .Fortran("subsparsefull",
                          nrow,dcheck(A@entries),A@colindices,
                          A@rowpointers,b=dcheck(B),PACKAGE = "spam"
                          )$b,nrow,ncol))
}

".spam.addsubsparsesparse" <-
function(A,B,s)
{
  nrow <- A@dimension[1]
  ncol <- A@dimension[2]
  if(ncol != B@dimension[2] || nrow != B@dimension[1])
    stop("non-conformable matrices")
  nzmax <- length(union(A@colindices+ncol*(rep(1:nrow,diff(A@rowpointers))-1),
                        B@colindices+ncol*(rep(1:nrow,diff(B@rowpointers))-1)))+1
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
                entries     = double(nzmax),
                colindices  = integer(nzmax),
                rowpointers = integer(nrow+1),
                as.integer(nzmax),
                ierr = integer(1),
                PACKAGE = "spam"
                )
  if(z$ierr != 0) stop("insufficient space for sparse matrix addition")
  nz <- z$rowpointers[nrow+1]-1
  return(new("spam",entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,dimension=c(nrow,ncol)))
}


".spam.elemul" <-
function(e1,e2)
{
  if(is.vector(e1)) {
    if(length(e1) == 1){
      if(e1==0) return(as.spam(0,nrow(e2),ncol(e2)))
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
      if(e2==0)   return(as.spam(0,nrow(e1),ncol(e1)))
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
                icheck(e1row),
                icheck(e1col),
                dcheck(e1@entries),
                icheck(e1@colindices),
                icheck(e1@rowpointers),
                dcheck(e2@entries),
                icheck(e2@colindices),
                icheck(e2@rowpointers),
                entries     = double(nnzmax),
                colindices  = integer(nnzmax),
                rowpointers = integer(e1row+1),
                integer(e1col),
                double(e1col),
                as.integer(nnzmax),
                ierr = integer(1),
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
	nzmax <- length(union(e1@colindices+ncol*(rep(1:nrow,diff(e1@rowpointers))-1),
                e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1)))+1
        z <- .Fortran("aedib",
                      icheck(nrow),
                      icheck(ncol),
                      icheck(1),
                      dcheck(e1@entries),
                      icheck(e1@colindices),
                      icheck(e1@rowpointers),
                      dcheck(e2@entries),
                      icheck(e2@colindices),
                      icheck(e2@rowpointers),
                      entries     = double(nzmax),
                      colindices  = integer(nzmax),
                      rowpointers = integer(nrow+1),
                      as.integer(nzmax),
                      integer(ncol),
                      double(ncol),
                      ierr = integer(1),
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
	nzmax <- length(union(e1@colindices+ncol*(rep(1:nrow,diff(e1@rowpointers))-1),
                e2@colindices+e2@dimension[2]*(rep(1:e2@dimension[1],diff(e2@rowpointers))-1)))+1
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
                entries     = double(nzmax),
                colindices  = integer(nzmax),
                rowpointers = integer(nrow+1),
                as.integer(nzmax),
		integer(ncol),
		double(ncol),
                ierr = integer(1),
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
	  function (x, i, j, drop=FALSE) 
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
  ic <- unique( z$ic[1:(z$nr+1)])
  return(new("spam",entries=z$c[1:nz],colindices=z$jc[1:nz],rowpointers=ic,
               dimension=as.integer(c(length(ic)-1,max(z$jc[1:nz])))))
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
      return( x[as.spam(rw)] )
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
                    allelem=double(nir),
                      PACKAGE="spam")$allelem)

  }
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
  if (all(diff(rw)==1) & all(diff(cl)==1)) {
    z <- .Fortran("submat",
                  nrow,
                  job=as.integer(1), # need values as well
                  i1=as.integer(rw[1]),
                  i2=as.integer(rw[length(rw)]),
                  j1=as.integer(cl[1]),
                  j2=as.integer(cl[length(cl)]),
                  dcheck(x@entries),icheck(x@colindices),icheck(x@rowpointers),
                  nr=as.integer(0),
                  nc=as.integer(0),
                  ao=dcheck(x@entries),jao=icheck(x@colindices),iao=icheck(x@rowpointers),
                  PACKAGE="spam")
    nz <- z$iao[z$nr+1]-1
    if (drop==TRUE && (z$nr==1 || z$nc==1))
      return(c(new("spam",entries=z$ao[1:nz],colindices=z$jao[1:nz],
                 rowpointers=z$iao[1:(z$nr+1)],dimension=c(z$nr,z$nc))))
    else
      return(new("spam",entries=z$ao[1:nz],colindices=z$jao[1:nz],
                 rowpointers=z$iao[1:(z$nr+1)],dimension=c(z$nr,z$nc)))
  }
  if (is.vector(rw) && is.vector(cl)) {
    nrw <- length(rw)   # length returns an integer, so is a product therof
    ncl <- length(cl)
    bnz <- nrw*ncl
    z <- .Fortran("getblock",
                  dcheck(x@entries),icheck(x@colindices),icheck(x@rowpointers),
                  nrw,as.integer(rw),
                  ncl,as.integer(cl),
                  bnz=bnz, b=double(bnz),jb=integer(bnz),ib=integer(nrw+1),
                  PACKAGE="spam")
    return(new("spam",entries=z$b[1:z$bnz],colindices=z$jb[1:z$bnz],
               rowpointers=z$ib[1:(nrw+1)],dimension=c(nrw,ncl)))
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
  nzmax <- as.integer(min(c(nrow*ncol), i@rowpointers[i@dimension[1]+1]+x@rowpointers[nrow+1]-2))
  if (length(value)!=  (i@rowpointers[i@dimension[1]+1]-1) )
    value <- rep(value, (i@rowpointers[i@dimension[1]+1]-1) %/%length(value))
#   cat(length(value))#@@#
  z <- .Fortran("subass",
                nrow,ncol,
                dcheck(x@entries),      x@colindices,    x@rowpointers,
                b=as.double(value),  bj=i@colindices, bi=i@rowpointers,
                c=double(nzmax),jc=integer(nzmax),ic=integer(nrow+1),
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
    rw <- rw[ord,]
    bia <- .Fortran("constructia",
                    nrow,
                    rowpointers=integer(nrow+1),
                    ir=as.integer(rw[,1]),
                    PACKAGE="spam")$rowpointers
    nzmax <- as.integer(min(nrow*ncol, nir+x@rowpointers[nrow+1]+1)+1)
    # new("spam",entries=value[ord],colindices=as.integer(rw[,2]),rowpointers=bia,c(nrow,ncol))
    z <- .Fortran("subass",
                  nrow,ncol,
                  dcheck(x@entries), x@colindices, x@rowpointers,
                  b=as.double(value[ord]),
                  bj=as.integer(rw[,2]),
                  bi=bia,
                  c=double(nzmax),jc=integer(nzmax),ic=integer(nrow+1),
                  nzmax=nzmax,
                  PACKAGE="spam")
    cnz <- z$ic[nrow+1]-1
    return(new("spam",entries=z$c[1:cnz],colindices=z$jc[1:cnz],
               rowpointers=z$ic[1:(nrow+1)],dimension=c(nrow,ncol)))
    
  }
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
    
    bia <- numeric(nrow)  # bia has size of nrow + 1
    bia[rw] <- ncl        # in each row we have ncl new objects
    bia <- as.integer(c(1,cumsum(bia)+1))
                
    # we construct now a sparse matrix containing the "value" at positions rw and cl.
    # then we use the subassign function.
    nzmax <- as.integer(min(c(nrow*ncol), bnz+x@rowpointers[nrow+1]-1)+2)
    # new("spam",entries=value,colindices=rep(sort(as.integer(cl)),nrw),rowpointers=bia,c(nrow,ncol))
    z <- .Fortran("subass",
                  nrow,ncol,
                  dcheck(x@entries), x@colindices ,x@rowpointers,
                  b=value,
                  bj=rep(sort(as.integer(cl)),nrw),
                  bi=bia,
                  c=double(nzmax),jc=integer(nzmax),ic=integer(nrow+1),
                  nzmax=nzmax,
                  PACKAGE="spam")
    cnz <- z$ic[nrow+1]-1
    return(new("spam",entries=z$c[1:cnz],colindices=z$jc[1:cnz],
               rowpointers=z$ic[1:(nrow+1)],dimension=c(nrow,ncol)))
  }
  stop("invalid or not-yet-implemented 'spam' subsetting")
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
                  y=double(nrow),
                  dcheck(y@entries),
                  y@colindices,
                  y@rowpointers,
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
                  y=double(nrow),
                  dcheck(x@entries),
                  x@colindices,
                  x@rowpointers,
                  PACKAGE = "spam")$y
    dim(z) <- c(nrow,1)
    return(z)
  }
  if (is.matrix(y)) y <- as.spam(y)
  else  if(is.matrix(x)) x <- as.spam(x)


  #matrix multiply two sparse spam matrices

  xn <- icheck(x@dimension[1])
  xm <- icheck(x@dimension[2])
  yl <- icheck(y@dimension[2])
  if(xm != icheck(y@dimension[1]))
    stop("matrices not conformable for multiplication")

  z <- .Fortran("amubdg",
                xn,xm,yl,
                x@colindices,x@rowpointers,
                y@colindices,y@rowpointers,
                integer(xn),
                nz = integer(1),
                integer(yl),
                PACKAGE = "spam")
  nzmax <- z$nz
  z <- .Fortran("amub",
                xn,yl,
                as.integer(1),
                dcheck(x@entries), x@colindices, x@rowpointers,
                dcheck(y@entries), y@colindices, y@rowpointers,
                entries = double(nzmax), colindices = integer(nzmax), rowpointers = integer(xn+1),
                as.integer(nzmax),
                integer(yl),
                ierr = integer(1),
                PACKAGE = "spam")
  nz <- z$rowpointers[xn+1]-1
  if(z$ierr != 0) stop("insufficient space for sparse matrix multiplication")
  
      
  if(nz==0){#trap zero matrix
    z$entries <- 0
    z$colindices <- as.integer(1)
    z$rowpointers <- as.integer(c(1,rep(2,xn)))
  }  else  z <- .Fortran("sortrows",
                         xn,entries=z$entries[1:nz],colindices=z$colindices[1:nz],rowpointers=z$rowpointers,
                         PACKAGE = "spam")
  
  return(new("spam",entries=z$entries,colindices=z$colindices,
             rowpointers=z$rowpointers,dimension=icheck(c(xn,yl))))
}



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

setMethod("Math","spam", function(x){ x@entries <- callGeneric(x@entries);x })
setMethod("Math2",signature(x = "spam", digits = "numeric"),
          function(x, digits){ x@entries <- callGeneric(x@entries, digits = digits);x })

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
setMethod("Summary","spam", function(x,...,na.rm=FALSE){ callGeneric(x@entries,...,na.rm=FALSE) })


setMethod("%*%",signature(x="spam",y="spam"),    .spam.matmul)
setMethod("%*%",signature(x="spam",y="matrix"),  .spam.matmul)
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
                      diag=double(nrow),
                      idiag=integer(nrow),
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
                  diag=double(nrow),
                  idiag=integer(nrow),
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
  switch(type,
         sup = max(abs(x)),
         l1 = sum(abs(x)),
         sqrt(sum(x^2))
         )
}

setGeneric("norm",function(x, type,...)standardGeneric("norm"))
setMethod("norm","spam", function(x, type = "sup", ...){
  switch(type,
         sup = max(abs(x@entries)),
         l1 = sum(abs(x@entries)),
         sqrt(sum(x@entries^2))
         )
}
          )

#setMethod("norm","numeric", function(x, type = "sup", ...){
#  switch(type,
#         sup = max(abs(x)),
#         HS = sqrt(sum(x^2)),
#         l1 = sum(abs(x))
#         )
#}
#          )

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
            length(breaks), code = integer(length(zvals)), (TRUE),
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
      points( xx[rep((1:nrow(z)),diff(z@rowpointers))], yy[z@colindices], pch='.', cex=cex*.Spam$cex/(ncol(z)+nrow(z)),
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
  if (nrow*ncol < .Spam$imagesize) {
    z <- numeric(nrow*ncol)
    dim(z) <- c(nrow,ncol)
    z[cbind(rep(nrow:1,diff(x@rowpointers)),x@colindices)] <- -1
    image.default(x=1:ncol,y=-(nrow:1),t(z),
                  axes=FALSE, col=col, xlab=xlab, ylab=ylab, ...) 
  } else {
    if (missing(cex)) {
      warning("default value for 'cex' in 'display' might not be the optimal choice", call.=FALSE)
      cex <- 1
    }
    plot( x@colindices, rep(-(1:nrow),diff(x@rowpointers)), pch='.', cex=cex*.Spam$cex/(ncol+nrow),
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
  tmp <- x[,1:2] # extract the first two columns
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
                icheck(nrow),
                entries=dcheck(y@entries),
                icheck(y@colindices),
                icheck(y@rowpointers),
                as.double(x),
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
                nrow = icheck(nrow),
                n = minrc,
                a = c(y@entries,double(minrc)),
                ja = c(y@colindices,integer(minrc)),
                ia = y@rowpointers,
                diag = as.double(x),
                iw = integer(nrow),    # just to be sure
                info = integer(nrow+1),
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


