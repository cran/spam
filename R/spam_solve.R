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

setClass("spam.chol.NgPeyton",representation(nrow="numeric",nnzlindx="numeric",
	nsuper="numeric",lindx="numeric",xlindx="numeric",nnzl="numeric",
	lnz="numeric",xlnz="numeric",invp="numeric",perm="numeric",
	xsuper="numeric"))


########################################################################

print.spam.chol.NgPeyton <- function(x,...) {
  cat(paste("(Upper) Cholesky factor of dimension ", x@nrow,
                "x", x@nrow, " with ",x@nnzl," nonzero elements.", sep = ""), fill=TRUE)
  cat("  (The object is supposed to be used with: 'as.spam', 'backsolve' or 'forwardsolve'.)\n", fill=TRUE)
  cat("Class 'spam.chol.NgPeyton'\n")
  invisible(x)
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
  invisible(object)
        })

setMethod("print","spam.chol.NgPeyton", print.spam.chol.NgPeyton)
setMethod("summary","spam.chol.NgPeyton",summary.spam.chol.NgPeyton)

auxiliarycholNgPeyton <- function(x,nrow,memory){
  nnzdmax <- x@rowpointers[nrow+1]-1
  nnzdsm <- nnzdmax + nrow + 1
  iwmax <- 7*nrow+3
  level <- 8

  if(is.null(memory$nsubmax))
    nsubmax <- nnzdmax
  else {
    nsubmax <- memory$nsubmax
    memory$nsubmax <- NULL
  }
  if(is.null(memory$nnzlmax))
    nnzlmax <- max(4*nnzdmax,floor(.2*nnzdmax^1.3))
  else {
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
  
  return(.Fortran("cholmod",
                nrow = icheck(nrow),
                nnzdmax = as.integer(nnzdmax),
                d =  dcheck(x@entries),
                jd = x@colindices,
                id = x@rowpointers,
                nnzdsm = as.integer(nnzdsm),
                dsub = double(nnzdsm),
                jdsub = integer(nnzdsm),
                nsub = integer(1),
                nsubmax = as.integer(nsubmax),
                lindx = integer(nsubmax),
                xlindx = integer(nrow+1),
                nsuper = integer(1),
                nnzlmax = as.integer(nnzlmax),
                lnz = double(nnzlmax),
                xlnz = integer(nrow+1),
                invp = integer(nrow),
                perm = integer(nrow),
                iwmax = as.integer(iwmax),
                iwork = integer(iwmax),
                colcnt = integer(nrow),
                snode = integer(nrow),
                xsuper = integer(nrow+1),
                split = integer(nrow),
                tmpmax = as.integer(tmpmax),
                tmpvec = double(tmpmax),
                cachsz = as.integer(cache),
                level = as.integer(level),
                ierr = integer(1),
                PACKAGE = "spam") )
}







"as.spam.chol.NgPeyton" <- function(x, eps = .Spam$eps){
  ja <- .Fortran('calcja',
                 x@nrow, x@nsuper, x@xsuper, x@lindx, x@xlindx, x@xlnz,
                 xja=as.integer(numeric(x@nnzl)),
                  PACKAGE = "spam")$xja
  return(new('spam',entries=x@lnz, colindices=ja, rowpointers=x@xlnz, dimension=c(x@nrow,x@nrow)))
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


"chol" <- function(x, ...) UseMethod("chol")
"chol.default" <- base::chol
setGeneric("chol")
setMethod("chol","matrix",base::chol)
setMethod("chol","spam", function(x, #pivot = FALSE,
        method="NgPeyton",
        ordering="MinDeg",                          
        memory=list(),
	eps = .Spam$eps, ...){

  nrow <- x@dimension[1]
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)
  if(norm(t(x)-x,type='sup') > 2*eps) stop("Input matrix to 'chol' not symmetric")

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

  z <- auxiliarycholNgPeyton(x,nrow,memory)
  
  if (z$ierr != 0){
    if(z$ierr == 9) mess <- "singularity problem"
    else if(z$ierr == 4) mess <- paste("Increase nnzlmax (currently set to ",z$nnzlmax,')',sep='')
    else if(z$ierr == 5) mess <- paste("Increase nsubmax (currently set to ",z$nsubmax,')',sep='')
    else if(z$ierr %in% c(8,10)) mess <- paste("Increase tmpmax (currently set to ",z$tmpmax,')',sep='')
    else mess <- "insufficient space"
    stop(mess)
  }
  nnzl <- z$xlnz[length(z$xlnz)]-1
  nnzlindx <- z$nsub

  
  invisible(new("spam.chol.NgPeyton",nrow=z$nrow,nnzlindx=nnzlindx,
      nsuper=z$nsuper,lindx=z$lindx[1:nnzlindx],xlindx=z$xlindx,
      nnzl=icheck(nnzl),lnz=z$lnz[1:nnzl],xlnz=z$xlnz,invp=z$invp,
      perm=z$perm,xsuper=z$xsuper))
})

setMethod("solve","spam",
function (a, b, ...) {
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
    z <- .Fortran("bckslv", m = icheck(nr), nnzlindx = icheck(a@nnzlindx),
                  icheck(a@nsuper), as.integer(p), icheck(a@lindx),
                  icheck(a@xlindx), icheck(a@nnzl), dcheck(a@lnz),
                  icheck(a@xlnz), icheck(a@invp), icheck(a@perm),
                  icheck(a@xsuper), double(nr), sol = double(nr*p), as.double(b),
                  PACKAGE = "spam")$sol
  else z <- backsolve(a,forwardsolve( t(a),b))
  
  if ( p!=1)    dim(z) <- c(nr,p)
  return( z)
})

setMethod("backsolve","spam.chol.NgPeyton",
function(r, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
#       r: chol.spam.NgPeyton structure as returned by chol.
#       x: rhs  may be a matrix in dense form
#  subroutine backsolve(m,nsuper,nrhs,lindx,xlindx,lnz,xlnz,xsuper,b)
        m <- r@nrow
        if(!is.matrix(x)) x <- as.matrix(x)
        if(nrow(x)!=m)stop("Cholesky factor not conformable with x")
        p <- ncol(x)
        if (!.Spam$bcksl) {
        z <- .Fortran("backsolve", m = icheck(m), 
                icheck(r@nsuper), as.integer(p),  icheck(r@lindx),
                icheck(r@xlindx), dcheck(r@lnz),  icheck(r@xlnz),  
                icheck(r@xsuper), sol =  as.double(x),
                      PACKAGE="spam")$sol
      }else{
        z <- .Fortran("bckslb", m = icheck(m), nnzlindx = icheck(r@nnzlindx),
                icheck(r@nsuper), as.integer(p),  icheck(r@lindx),
                icheck(r@xlindx), icheck(r@nnzl), dcheck(r@lnz),
                icheck(r@xlnz),   icheck(r@invp), icheck(r@perm),
                icheck(r@xsuper), double(m), sol = double(m*p), as.double(x),
                      PACKAGE="spam")$sol
      }
        if (p>1)
          dim(z) <- c(nrow(x),p)
        return(z)
        })

setMethod("forwardsolve","spam.chol.NgPeyton",
function(l, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
#       l: chol.spam.NgPeyton structure as returned by chol
#       x: rhs  may be a matrix in dense form
        m <- l@nrow
        if(!is.matrix(x)) x <- as.matrix(x)
        if(nrow(x)!=m) stop("Cholesky factor not conformable with x")
        p <- ncol(x)
        if (!.Spam$bcksl) {
        z <- .Fortran("forwardsolve", m = icheck(m), 
                icheck(l@nsuper), as.integer(p), icheck(l@lindx),
                icheck(l@xlindx), dcheck(l@lnz), icheck(l@xlnz), 
                icheck(l@xsuper), sol = as.double(x),
                      PACKAGE="spam")$sol
      }else{
        z <- .Fortran("bckslf", m = icheck(m), nnzlindx = icheck(l@nnzlindx),
                icheck(l@nsuper), as.integer(p),  icheck(l@lindx),
                icheck(l@xlindx), icheck(l@nnzl), dcheck(l@lnz),
                icheck(l@xlnz),   icheck(l@invp), icheck(l@perm),
                icheck(l@xsuper), double(m), sol = double(m*p), as.double(x),
                      PACKAGE="spam")$sol
      }
        if (p>1)
          dim(z) <- c(nrow(x),p)
        return(z)
        })

######################################################################
######################################################################

determinant.spam <- function(x,  logarithm = TRUE,
        method="NgPeyton",
        ordering="MinDeg",                          
        memory=list(),
	eps = .Spam$eps, ...){
  
  nrow <- x@dimension[1]
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)
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

  z <- auxiliarycholNgPeyton(x,nrow,memory)

 
  if (z$ierr != 0){
    if(z$ierr == 9) {
      warning("singularity problem or matrix not positive definite",call.=FALSE)
      logdet$modulus <- NA
    }
    else if(z$ierr == 4) stop("Increase nnzlmax with 'NgPeyton' method",call.=FALSE)
    else if(z$ierr == 5) stop("Increase nsubmax with 'NgPeyton' method",call.=FALSE)
    else if(z$ierr %in% c(8,10)) stop("Increase tmpmax with 'NgPeyton' method",call.=FALSE)
    else stop("Insufficient space for 'NgPeyton' method",call.=FALSE)
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
