########################################################################
########################################################################
#
# Contains routines linked to solving linear systems. Namely:
#    chol, solve, backsolve, forwardsolve
#    determinant                (the later because it is based on chol)
#
#    As well as associated S4 elements.
# 
########################################################################
########################################################################

setClass("spam.chol.NgPeyton",
         representation(nrow="integer",nnzlindx="integer",
                        nsuper="integer",lindx="integer",xlindx="integer",nnzl="integer",
                        lnz="numeric",xlnz="integer",invp="integer",perm="integer",
                        xsuper="integer"),
         # the prototype corresponds to the cholesky of '1'
         prototype=prototype(nrow=as.integer(1),nnzlindx=as.integer(1),
           nsuper=as.integer(1),lindx=as.integer(1),xlindx=as.integer(c(1,2)),
           nnzl=as.integer(1),lnz=1.0,xlnz=as.integer(c(1,2)),
           invp=as.integer(1),perm=as.integer(1),xsuper=as.integer(c(1,2))
           )
         )


########################################################################

print.spam.chol.NgPeyton <- function(x,...) {
  cat(paste("(Upper) Cholesky factor of dimension ", x@nrow,
                "x", x@nrow, " with ",x@nnzl," nonzero elements.", sep = ""), fill=TRUE)
  cat("  (The object is supposed to be used with: 'as.spam', 'backsolve' or 'forwardsolve'.)\n", fill=TRUE)
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(NULL)
}
summary.spam.chol.NgPeyton <- function(object,...) {
  dens <- object@nnzl/(object@nrow^2) * 100
  cat(paste("(Upper) Cholesky factor of class 'spam.chol.NgPeyton' of dimension ", object@nrow,
                "x", object@nrow, " with ",object@nnzl," nonzero elements.", sep = ""), fill=TRUE)
  cat(paste("  Density of the matrix is ", signif(dens, 3),
        "% (nz=", object@nnzl, ").\n", sep = ""))
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(object)
}

setMethod("show","spam.chol.NgPeyton", function(object) {
  cat(paste("(Upper) Cholesky factor of dimension ", object@nrow,
                "x", object@nrow, " with ",object@nnzl," nonzero elements.", sep = ""), fill=TRUE)
  cat("  (The object is supposed to be used with: 'as.spam', 'backsolve' or 'forwardsolve'.)\n", fill=TRUE)
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(NULL)
        })

setMethod("print","spam.chol.NgPeyton", print.spam.chol.NgPeyton)
setMethod("summary","spam.chol.NgPeyton",summary.spam.chol.NgPeyton)



"as.spam.chol.NgPeyton" <- function(x, eps = .Spam$eps){
  newx <- new("spam")
  slot(newx,"entries",check=FALSE) <- x@lnz
  slot(newx,"colindices",check=FALSE) <- .Fortran('calcja',
                           x@nrow, x@nsuper, x@xsuper, x@lindx, x@xlindx, x@xlnz,
                           xja=vector("integer",x@nnzl),
                           NAOK = !.Spam$safemode,
                           DUP=FALSE,
                           PACKAGE = "spam")$xja
  slot(newx,"rowpointers",check=FALSE) <- x@xlnz
  slot(newx,"dimension",check=FALSE) <- c(x@nrow,x@nrow)
  return(newx)
}

setMethod("as.spam","spam.chol.NgPeyton", as.spam.chol.NgPeyton)


"backsolve" <- function(r,x, ...) UseMethod("backsolve")
#"backsolve.default" <- base::backsolve
setGeneric("backsolve")
setMethod("backsolve","matrix",base::backsolve)

"forwardsolve" <- function(l,x, ...) UseMethod("forwardsolve")
#"forwardsolve.default" <- base::forwardsolve
setGeneric("forwardsolve")
setMethod("forwardsolve","matrix",base::forwardsolve)

"ordering.default" <- function(x,inv=FALSE) stop('Operation not defined form this class')

#ordering <- function(x,...) stop('Operation not defined form this class')
#setGeneric("ordering")
setGeneric("ordering",function(x,inv=FALSE)standardGeneric("ordering"))

setMethod("ordering","spam.chol.NgPeyton",function(x,inv=FALSE)
          {
            if (inv) return(x@invp) else return(x@perm) })



setMethod("ordering","matrix",function(x,inv=FALSE)
          {
            if (dim(x)[1]!=dim(x)[2])
              stop("ordering is defined for square matrices only")
            if(inv)return(dim(x)[1]:1) else return(1:dim(x)[1]) })

setMethod("ordering","spam",function(x,inv=FALSE)
          {
            if (dim(x)[1]!=dim(x)[2])
              stop("ordering is defined for square matrices only")
            if(inv)return(dim(x)[1]:1) else return(1:dim(x)[1]) })


if(paste(R.version$major, R.version$minor, sep=".") < "2.6") 
{
 
"chol" <- function(x, ...) UseMethod("chol")
"chol.default" <- base::chol
setGeneric("chol")
setMethod("chol","matrix",base::chol)

}

chol.spam <- function(x, #pivot = FALSE,
                                  method="NgPeyton",
                                  ordering="MinDeg",                          
                                  memory=list(),
                                  eps = .Spam$eps, ...){

  nrow <- x@dimension[1]
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)

  if(.Spam$safemode) {
    if(norm(t(x)-x,type='sup') > (2+eps)*eps) stop("Input matrix to 'chol' not symmetric (up to (2+eps)*eps in 'sup'-norm")
  }
  
  if (method != "NgPeyton")
    warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'",
                     method), domain = NA)

  if (is.numeric(ordering)) {
    if (length(ordering) != nrow)
      stop("Permutation defined in 'ordering' is of wrong length")
    warning("manual ordering is not supported. Using 'MinDeg'", domain = NA)
  }
  if (ordering != "MinDeg")
    warning(gettextf("ordering = '%s' is not supported. Using 'MinDeg'",
                     method), domain = NA)

#  if (!is.null(ordering)) 

  nnzdmax <- x@rowpointers[nrow+1]-1
  nnzdsm <- nnzdmax + nrow + 1
  iwmax <- 7*nrow+3
  level <- 8

  if(is.null(memory$nsubmax))    nsubmax <- nnzdmax  else {
    nsubmax <- max(memory$nsubmax,nnzdmax)
    memory$nsubmax <- NULL
  }
  if(is.null(memory$nnzlmax))    nnzlmax <- max(4*nnzdmax,floor(.2*nnzdmax^1.3))  else {
    nnzlmax <- memory$nnzlmax
    memory$nnzlmax <- NULL
  }
  
  if(is.null(memory$tmpmax))    tmpmax <- 500*nrow  else {
    tmpmax <- memory$tmpmax 
    memory$tmpmax <- NULL
  }
  if(is.null(memory$cache))    cache <- 64  else {
    cache <- memory$cache 
    memory$cache <- NULL
  }

  if (length( memory)>0 )
    warning("The component(s) ", paste("'",names(memory),"'",sep='',collapse=",")," of the argument 'memory'\npassed to function 'chol' are not meaningful and hence ignored.",call.=FALSE)
  
  z <- .Fortran("cholmod",
                   nrow = nrow,                   #1
                   nnzdmax = as.integer(nnzdmax),
                   d =  dcheck(x@entries),
                   jd = x@colindices,
                   id = x@rowpointers,
                   nnzdsm = as.integer(nnzdsm),
                   dsub = vector("double",nnzdsm),
                   jdsub = vector("integer",nnzdsm),
                   nnzlindx = vector("integer",1),            #9 formerly nnzlindx <- z$nsub
                   nsubmax = as.integer(nsubmax),
                   lindx = vector("integer",nsubmax),     #11
                   xlindx = vector("integer",nrow+1),     #
                   nsuper = vector("integer",1),          #
                   nnzlmax = as.integer(nnzlmax),#
                   lnz = vector("double",nnzlmax),        #
                   xlnz = vector("integer",nrow+1),     #
                   invp = vector("integer",nrow),       #
                   perm = vector("integer",nrow),       #18
                   iwmax = as.integer(iwmax),
                   iwork = vector("integer",iwmax),
                   colcnt = vector("integer",nrow),
                   snode = vector("integer",nrow),
                   xsuper = vector("integer",nrow+1),   #23
                   split = vector("integer",nrow),
                   tmpmax = as.integer(tmpmax),
                   tmpvec = vector("double",tmpmax),
                   cachsz = as.integer(cache),
                   level = as.integer(level),
                   ierr = vector("integer",1),           #29
                NAOK = .Spam$safemode,
                   DUP=!FALSE,PACKAGE = "spam")
  
  z[c(2:8,10,19:22,24:28)] <- NULL
  
  if(z$ierr == 4) stop(paste("Increase 'nnzlmax' with 'NgPeyton' method (currently set to ",nnzlmax,")",sep=""),call.=FALSE)
  if(z$ierr == 5) stop(paste("Increase 'nsubmax' with 'NgPeyton' method (currently set to ",nsubmax,")",sep=""),call.=FALSE)
  if(z$ierr %in% c(8,10)) stop(paste("Increase 'tmpmax' with 'NgPeyton' method (currently set to ",tmpmax,")",sep=""),call.=FALSE)
  if(z$ierr>0&z$ierr!=9) stop("Insufficient space for 'NgPeyton' method",call.=FALSE)
  if(z$ierr == 9) stop("singularity problem") 
  nnzl <- as.integer(z$xlnz[length(z$xlnz)]-1)

  newx <- new("spam.chol.NgPeyton")
  slot(newx,"nrow",check=FALSE) <- nrow
  slot(newx,"nnzlindx",check=FALSE) <- z$nnzlindx
  slot(newx,"nsuper",check=FALSE) <- z$nsuper
  slot(newx,"lindx",check=FALSE) <- z$lindx[1:z$nnzlindx]
  slot(newx,"xlindx",check=FALSE) <- z$xlindx
  slot(newx,"nnzl",check=FALSE) <- nnzl
  slot(newx,"lnz",check=FALSE) <- z$lnz[1:nnzl]
  slot(newx,"xlnz",check=FALSE) <- z$xlnz
  slot(newx,"invp",check=FALSE) <- z$invp
  slot(newx,"perm",check=FALSE) <- z$perm
  slot(newx,"xsuper",check=FALSE) <- z$xsuper
  
  
  invisible(newx)
}

solve.spam <- function (a, b, ...) {
  nr <- a@dimension[1]
  nc <- a@dimension[2]
  if (nc != nr)      stop("'A' must be square")
  a <- chol(a,...)

  if (missing(b)) {
    b <- diag(1, nc)
    missing.b <- TRUE
  }
  else {
    missing.b <- FALSE
    if(!is.matrix(b)) b <- as.matrix(b)
  }
  p <- ncol(b)
  if (is(a,"spam.chol.NgPeyton"))
      # The following is a fast way to perform:
      #     z <- backsolve(a,forwardsolve( t(a),b))
    z <- .Fortran("bckslv", m = nr, nnzlindx = a@nnzlindx,
                  a@nsuper, as.integer(p), a@lindx,
                  a@xlindx, a@nnzl, dcheck(a@lnz),
                  a@xlnz, a@invp, a@perm,
                  a@xsuper, double(nr), sol = vector("double",nr*p), as.double(b),DUP=FALSE,
                  NAOK = !.Spam$safemode,
                  PACKAGE = "spam")$sol
  else z <- backsolve(a,forwardsolve( t(a),b))
  
  if ( p!=1)    dim(z) <- c(nr,p)
  return( z)
}
backsolve.spam <- function(r, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
#       r: spam.chol.NgPeyton structure as returned by chol.
#       x: rhs  may be a matrix in dense form
#  subroutine backsolve(m,nsuper,nrhs,lindx,xlindx,lnz,xlnz,xsuper,b)
  m <- r@nrow
  if(is.vector(x)) {
    n <- length(x)
    p <- as.integer(1)
  } else {
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
  }
  if(n!=m)stop("Cholesky factor not conformable with x")
  if (!.Spam$bcksl) {
    z <- .Fortran("backsolve", m = m, 
                  r@nsuper, as.integer(p),  r@lindx,
                  r@xlindx, dcheck(r@lnz),  r@xlnz,  
                  r@xsuper, sol =  as.double(x),#DUP=FALSE,
                  NAOK = !.Spam$safemode,
                  PACKAGE="spam")$sol
  }else{
    z <- .Fortran("bckslb", m = m, nnzlindx = icheck(r@nnzlindx),
                  r@nsuper, as.integer(p),  r@lindx,
                  r@xlindx, r@nnzl, dcheck(r@lnz),
                  r@xlnz,   r@invp, r@perm,
                  r@xsuper, double(m), sol = vector("double",m*p), as.double(x),DUP=FALSE,
                  NAOK = !.Spam$safemode,
                  PACKAGE="spam")$sol
  }
  if (p>1)
    dim(z) <- c(m,p)
  return(z)
}
forwardsolve.spam <- function(l, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
#       l: spam.chol.NgPeyton structure as returned by chol
#       x: rhs  may be a matrix in dense form
  m <- l@nrow
  if(is.vector(x)) {
    n <- length(x)
    p <- as.integer(1)
  } else {
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
  }
  if(n!=m) stop("Cholesky factor not conformable with x")
  if (!.Spam$bcksl) {
    z <- .Fortran("forwardsolv", m = m, 
                  l@nsuper, as.integer(p), l@lindx,
                  l@xlindx, dcheck(l@lnz), l@xlnz, 
                  l@xsuper, sol = as.double(x),#DUP=FALSE,
                  NAOK = !.Spam$safemode,
                  PACKAGE="spam")$sol
  }else{
    z <- .Fortran("bckslf", m = m, nnzlindx = l@nnzlindx,
                  l@nsuper, as.integer(p),  l@lindx,
                  l@xlindx, l@nnzl, dcheck(l@lnz),
                  l@xlnz,   l@invp, l@perm,
                  l@xsuper, double(m), sol = vector("double",m*p), as.double(x),DUP=FALSE,
                  NAOK = !.Spam$safemode,
                  PACKAGE="spam")$sol
  }
  if (p>1)
    dim(z) <- c(m,p)
  return(z)
}



setMethod("chol","spam", chol.spam)
setMethod("solve","spam",solve.spam)

setMethod("backsolve","spam.chol.NgPeyton",backsolve.spam)
setMethod("forwardsolve","spam.chol.NgPeyton",forwardsolve.spam)

######################################################################
######################################################################

determinant.spam <- function(x, logarithm = TRUE, method="NgPeyton",
                             ordering="MinDeg", memory=list(),eps =
  .Spam$eps, ...){
  
  nrow <- x@dimension[1]
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)

  if (.Spam$safemode)
    if(norm(t(x)-x,type='sup') > 2*eps) stop("Input matrix to chol() not symmetric",call.=FALSE)
  
  if (method != "NgPeyton")
    warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'",
                     method), domain = NA)
  
  if (is.numeric(ordering)) {
    if (length(ordering) != nrow)
      stop("Permutation defined in 'ordering' is of wrong length",call.=FALSE)
    warning("manual ordering is not supported. Using 'MinDeg'", domain = NA)
  }
  if (ordering != "MinDeg")
    warning(gettextf("ordering = '%s' is not supported. Using 'MinDeg'",
                     method), domain = NA)
  


  logdet <- list()

  nnzdmax <- x@rowpointers[nrow+1]-1
  nnzdsm <- nnzdmax + nrow + 1
  iwmax <- 7*nrow+3
  level <- 8

  if(is.null(memory$nsubmax))  nsubmax <- nnzdmax  else {
    nsubmax <- max(memory$nsubmax,nnzdmax)
    memory$nsubmax <- NULL
  }
  if(is.null(memory$nnzlmax))    nnzlmax <- max(4*nnzdmax,floor(.2*nnzdmax^1.3))  else {
    nnzlmax <- memory$nnzlmax
    memory$nnzlmax <- NULL
  }
  
  if(is.null(memory$tmpmax))
    tmpmax <- 500*nrow
  else {
    tmpmax <- memory$tmpmax 
    memory$tmpmax <- NULL
  }
  if(is.null(memory$cache))
    cache <- 64
  else {
    cache <- memory$cache 
    memory$cache <- NULL
  }

  if (length( memory)>0 )
    warning("The components ", paste("'",names(memory),"'",sep='',collapse=",")," of the argument 'memory'\npassed to 'chol' are not meaningful and hence ignored.",call.=FALSE)
  
  z <- .Fortran("cholmod",
                   nrow = nrow,                   #1
                   nnzdmax = as.integer( nnzdmax),
                   d =  dcheck( x@entries),
                   jd = x@colindices,
                   id = x@rowpointers,
                   nnzdsm = as.integer( nnzdsm),
                   dsub = vector("double",nnzdsm),
                   jdsub = vector("integer",nnzdsm),
                   nnzlindx = vector("integer",1),            #9 formerly nnzlindx <- z$nsub
                   nsubmax = as.integer( nsubmax),
                   lindx = vector("integer",nsubmax),     #11
                   xlindx = vector("integer",nrow+1),     #
                   nsuper = vector("integer",1),          #
                   nnzlmax = as.integer( nnzlmax),#
                   lnz = vector("double",nnzlmax),        #
                   xlnz = vector("integer",nrow+1),     #
                   invp = vector("integer",nrow),       #
                   perm = vector("integer",nrow),       #18
                   iwmax = as.integer( iwmax),
                   iwork = vector("integer",iwmax),
                   colcnt = vector("integer",nrow),
                   snode = vector("integer",nrow),
                   xsuper = vector("integer",nrow+1),   #23
                   split = vector("integer",nrow),
                   tmpmax = as.integer( tmpmax),
                   tmpvec = vector("double",tmpmax),
                   cachsz = as.integer( cache),
                   level = as.integer( level),
                   ierr = vector("integer",1),           #29
                NAOK = !.Spam$safemode,
                   DUP=FALSE,PACKAGE = "spam")
  
  z[c(2:8,10,19:22,24:28)] <- NULL
  
  if(z$ierr == 4) stop(paste("Increase 'nnzlmax' with 'NgPeyton' method (currently set to ",nnzlmax,")",sep=""),call.=FALSE)
  if(z$ierr == 5) stop(paste("Increase 'nsubmax' with 'NgPeyton' method (currently set to ",nsubmax,")",sep=""),call.=FALSE)
  if(z$ierr %in% c(8,10)) stop(paste("Increase 'tmpmax' with 'NgPeyton' method (currently set to ",tmpmax,")",sep=""),call.=FALSE)
  if(z$ierr>0&z$ierr!=9) stop("Insufficient space for 'NgPeyton' method",call.=FALSE)
  if(z$ierr == 9) {
                                        # all other errors trapped 
      warning("singularity problem or matrix not positive definite",call.=FALSE)
      logdet$modulus <- NA
   } else{
    k <- z$xlnz
    tmp <- 2* sum( log(z$lnz[k[-length(k)]]))
    if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)
  }

  attr(logdet$modulus,"logarithm") <- logarithm
  
  logdet$sign <- ifelse(z$ierr == 9,NA,1)
  attr(logdet,"class") <- "det"
  
  return(logdet)
}

determinant.spam.chol.NgPeyton <- function(x, logarithm = TRUE,...)
{
  logdet <- list()

  k <- x@xlnz
  tmp <- sum( log(x@lnz[k[-length(k)]]))
  if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)
 
  attr(logdet$modulus,"logarithm") <- logarithm
  
  logdet$sign <- 1
  attr(logdet,"class") <- "det"
  
  return(logdet)
}


setMethod("determinant","spam", determinant.spam)
setMethod("determinant","spam.chol.NgPeyton", determinant.spam.chol.NgPeyton)

######################################################################
