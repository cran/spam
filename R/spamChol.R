# HEADER ####################################################
# This is file spam/R/spam_solve.R.                         #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



########################################################################
########################################################################
#
# Contains routines linked to solving spd linear systems. Namely:
#    chol.spam and chol.update
# Associated S4 elements and determinant are in different files.
#
#
########################################################################
########################################################################



# lindx=  colindices       of length nsub   = nnzR
# xlindx= colpointers                nsuper = nnzcolindices
# xlnz=   rowpointers
# snode=  snmember
# xsuper= supernodes
# c(... nnztmp,cachesize)= memory


#### Summary of memory usage in cholstepwise:
# Fortran      R
# m            nrow,ncol
# nsuper                     nnzcolindices
# nnzlmax                    nnzR
#
# Internally, we additionally assign several arrays (of length):
#  m: adj, colcnt
#  7*m+3: iwork
#  iwork: but 7*m+3 is used!



# TODO below, check with determinant.spam for use of vector_ds and potential switches
# from "rw" to "w" signature.
chol.spam <- function(x, pivot = "MMD",
                      method="NgPeyton",
                      memory=list(),
                      eps = getOption("spam.eps"), Rstruct=NULL, ..., verbose=FALSE){

  if (verbose) timeused <- proc.time()
  if (is(Rstruct,"spam.chol.NgPeyton")) invisible( update.spam.chol.NgPeyton(Rstruct,x,...))
  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)

  nrow <- x@dimension[1]
  nnzA <- x@rowpointers[nrow+1]-1L
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)

  if (any( diag.of.spam(x, nrow, nrow) < getOption("spam.eps")))
    stop("Input matrix to 'chol' not positive definite (up to eps)",call.=FALSE)
# base rule:
#  `nrow(x) * .Machine$double.neg.eps * max(diag(x)`.

    if (length(pivot)==1) {
        if (pivot==FALSE) {
            doperm <- 0
         } else if(pivot==TRUE) {
            doperm <- 1
        } else {
            doperm <- as.integer( switch(match.arg(pivot,c("MMD","RCM")),MMD=1,RCM=2))
        }
    } else  if (length(pivot)==nrow) {
        doperm <- 0
        if (getOption("spam.cholpivotcheck")) {
            checkpivot(pivot,nrow)
        }
    } else stop("'pivot' should be 'MMD', 'RCM' or a valid permutation")

  ### IMPROVEME get better parameter values
  nnzcfact <- c(5,1,5)
  nnzRfact <- c(5,1,2)
  # nnzcolindices = length of array holding the colindices
  if(is.null(memory$nnzcolindices))  {
   nnzcolindices <- ifelse((nnzA/nrow < 5), # very sparse matrix
                           max(1000,nnzA*(1.05*nnzA/nrow-3.8)),
                           nnzA)*nnzcfact[doperm+1]
  }else {
    nnzcolindices <- max(memory$nnzcolindices,nnzA)
    memory$nnzcolindices <- NULL
  }
  # nnzR = length of array holding the nonzero values of the factor
  if(is.null(memory$nnzR)) {
    nnzR <- min(max(4*nnzA,floor(.2*nnzA^1.3))*nnzRfact[doperm+1],nrow*(nrow+1)/2)
  } else {
    nnzR <- memory$nnzR
    memory$nnzR <- NULL
  }

#!4#
  if(is.null(memory$cache))     cache <- 512  else {
    cache <- memory$cache
    memory$cache <- NULL
  }

  if (length( memory)>0 )
    warning("The component(s) ", paste("'",names(memory),"'",sep='',collapse=","),
            " of the argument 'memory'\npassed to function 'chol' not meaningful and hence ignored.",call.=FALSE)


  if(getOption("spam.cholsymmetrycheck")) {
    test <- isSymmetric(x, tol = eps*100)
#  from help of isSymmetric:
#     isSymmetric(object, tol = 100 * .Machine$double.eps, ...)
    if (!isTRUE(test))
      stop("Input matrix to 'chol' not symmetric (up to 100*eps)",call.=FALSE)
  }

  if (method != "NgPeyton")
    warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'",
                     method), domain = NA)

    if (verbose) {
        cat("Factorizing matrix of dimension ", nrow, " with ", length(x@entries)," entries (",
            printSize(x),  "). \nReserved memory for factor with ", nnzR, " entries (",
            printSize(nnzR),").", sep='')
        }

    force64 <- getOption("spam.force64")
    # we elaborate on the sizes and check for 2^31-2 to allow one additional iteration per counter
    if(force64 || 7*nrow+3  > 2147483646 || nnzcolindices  > 2147483646 || nnzR  > 2147483646) {
        SS <- .format64()
    } else {
        SS <- .format32
    }
  #!3#


    if (length(pivot)==1) {
        if (pivot==FALSE) {
            pivot <- as.vector( seq_len(nrow), SS$type)
        } else if(pivot==TRUE) {
            pivot <- vector(SS$type,nrow)
        } else {
           pivot <- vector(SS$type, nrow)
        }
    } else  pivot <- as.vector( pivot, SS$type)



    if (verbose) {
        cat("\nWorking with",SS$name,"spam matrix")
        if (SS$type != class(x@dimension)) cat(", casting original matrix")

        ## note that object.size(x)  is not required as we have reading mode
        m1 <- object.size(x) +
            6*object.size(pivot)+
            (nnzcolindices +     # length lindx
             2*nrow +nnzA + 7*nrow+3)*   # work arrays {colcnt, split}, {adj, adjncy}  iwork
            ifelse(SS$signature=='double',8,4) +      # depending on 32/64-bit
            nnzR*8    # final entries
        cat(".\nApproximate amount of memory needed: ", printSize( m1),".", sep='')
                                        # gestimate for total use is:
#        m2 <- (8+6)*object.size(pivot)+
#            nnzcolindices*ifelse(SS$signature=='double',8,4) + nnzR*8
    }
  cholstep.intent <- c("r", "r", "r", "r",
                    "r", "r", "rw", "rw",
                    "rw", "r", "w", "w",    # 'w' for lindx xlindx  # _FOR SURE_
                    "rw", "rw", "w", "w",   # 'w' for lnz  xlnz
                    "rw", "rw", "r", "rw")

  cholstep.signature <- rep(SS$signature, 20)
  cholstep.signature[c(3,15)] <- "double"
  z <- .C64("cholstepwise",                # 11 huge objects...
            SIGNATURE=cholstep.signature,

            nrow = as.vector(nrow, mode=SS$type),
            nnzA = as.vector(nnzA, SS$type),
            d =  x@entries,                          #Size _D_nrow
            jd = as.vector(x@colindices, SS$type),   #Size nrow

            id = as.vector(x@rowpointers, SS$type),#5#Size nrow
            doperm = as.vector(doperm,    SS$type),
            invp = vector(SS$type,nrow),             #Size nrow
            perm = pivot,                            #Size nrow

            nnzlindx = vector(SS$type, 1), #9
            nnzcolindices = as.vector(nnzcolindices, SS$type),
            lindx = vector_dc(SS$type, nnzcolindices), #Size nnzcolindices
            xlindx = vector_dc(SS$type, nrow+1),       #Size nrow

            nsuper = vector(SS$type, 1), #13
            nnzR = as.vector(nnzR, SS$type),
            lnz = vector_dc( "double", nnzR),      #Size _D_nnzR
            xlnz = vector_dc( SS$type, nrow+1),    #Size nrow

            snode = vector(SS$type, nrow), #17     #Size nrow
            xsuper = vector(SS$type, nrow+1),      #Size nrow
            cachesize = as.vector(cache, SS$type),
            ierr = vector(SS$type, 1),

            INTENT=cholstep.intent,
            NAOK = getOption("spam.NAOK"),
            PACKAGE = SS$package)

  if (verbose) cat(".\n(Final size of R-Fortran transfered object: ", printSize( z),").", sep='')

  if(z$ierr == 1) stop("Singularity problem when calculating the Cholesky factor.")
  if(z$ierr == 6) stop("Inconsitency in the input",call.=FALSE)

  while( z$ierr>1) {
    if(z$ierr == 4) {
      tmp <- ceiling(nnzR*getOption("spam.cholincreasefactor")[1])
      warning("Increased 'nnzR' with 'NgPeyton' method\n",
                    "(currently set to ",tmp," from ",nnzR,")",call.=FALSE)
      nnzR <- tmp
    }
    if(z$ierr == 5) {
      tmp <- ceiling(nnzcolindices*getOption("spam.cholincreasefactor")[2])
      warning("Increased 'nnzcolindices' with 'NgPeyton' method\n",
         "(currently set to ",tmp," from ",nnzcolindices,")",call.=FALSE)
      nnzcolindices <- tmp
    }
    if(force64 || nnzcolindices  > 2147483647) {
        SS <- .format64()
    } else {
        SS <- .format32
    }
    cholstep.signature <- rep(SS$signature, 20)
    cholstep.signature[c(3,15)] <- "double"


    z <- .C64("cholstepwise",
              SIGNATURE=cholstep.signature,

              nrow = nrow,
              nnzA = x@rowpointers[nrow+1]-1,
              d =  x@entries,
              jd = x@colindices,

              id = x@rowpointers,
              doperm = doperm,
              invp = vector(SS$type,nrow),
              perm = pivot,

              nnzlindx = vector(SS$type,1),   # nnzR !!!
              nnzcolindices = nnzcolindices,
              lindx = vector(SS$type, nnzcolindices),
              xlindx = vector(SS$type, nrow+1),     #

              nsuper = 1,          #  number of supernodes
              nnzR = nnzR,#
              lnz = vector("double",nnzR),        #
              xlnz = vector(SS$type, nrow+1),     #

              snode = vector(SS$type, nrow),
              xsuper = vector(SS$type, nrow+1),
              cachesize = cache,
              ierr = 0L,

              INTENT=cholstep.intent,
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)

    if(z$ierr == 1) stop("Singularity problem when calculating the Cholesky factor.")
  }

  nnzR <- z$xlnz[length(z$xlnz)]-1

  newx <- new("spam.chol.NgPeyton")
  slot(newx,"entries",check=FALSE) <- z$lnz[1:nnzR]
  slot(newx,"colindices",check=FALSE) <- z$lindx[1:z$nnzlindx]
  slot(newx,"colpointers",check=FALSE) <- z$xlindx[1:(z$nsuper+1)]
  slot(newx,"rowpointers",check=FALSE) <- z$xlnz
  slot(newx,"dimension",check=FALSE) <- c(nrow,nrow)
  slot(newx,"pivot",check=FALSE) <- z$perm
  slot(newx,"invpivot",check=FALSE) <- z$invp
  slot(newx,"supernodes",check=FALSE) <- z$xsuper[1:(z$nsuper+1)]
  slot(newx,"snmember",check=FALSE) <- z$snode
  slot(newx,"memory",check=FALSE) <- c(nnzcolindices,z$nnzR,cache)
  slot(newx,"nnzA",check=FALSE) <- nnzA

    if (verbose) {
        cat("Final size of Cholesky factor (of class 'spam.chol.NgPeyton'): ")
        print( object.size(newx), units=c("Kb","Mb","Gb")[floor( log(as.numeric(m1),1024))])

        cat("Finished factorization (total time used:", (proc.time()-timeused)[1],"s).\n\n")
        }
  invisible(newx)
}



########################################################################


update.spam.chol.NgPeyton <- function(object,x,...){
    ## print("update.spam.chol.NgPeyton")
  nrow <- object@dimension[1]
  if (!is.spam(x))
    stop("Covariance should be a 'spam' object.")
  if ((x@rowpointers[nrow+1]-1) != object@nnzA)
    stop("Updated covariance entries do not match length of original one.")

    if( getOption("spam.force64") || .format.spam(object)$package != "spam" || .format.spam(x)$package != "spam" )
        SS <- .format64()
    else
        SS <- .format32

    u <- .C64("updatefactor",
     ##          subroutine updatefactor( m,nnzd,
     ## &     d,jd,id, invp,perm,
     ## &                lindx,xlindx, nsuper,lnz,xlnz,
     ## &                snode, xsuper,
     ## &                cachesize,ierr)
              SIGNATURE = c(SS$signature, SS$signature,
                            "double", SS$signature, SS$signature,
                            SS$signature, SS$signature, SS$signature, SS$signature, SS$signature,
                            "double", SS$signature, SS$signature, SS$signature, SS$signature, SS$signature ),
              nrow,
              object@nnzA,

              d =  x@entries,
              jd = x@colindices,
              id = x@rowpointers,

              object@invpivot,
              object@pivot,
              lindx = object@colindices,
              xlindx = object@colpointers,
              nsuper = length(object@supernodes)-1,

              entries = vector_dc( "double", length(object@entries)), #lnz
              rowpointers = object@rowpointers,#xlnz
              snode = object@snmember,
              xsuper = object@supernodes,
              cachesize = object@memory[3],
              ierr = 0,

              ## INTENT = c("r", "r",
              ##            "r", "r", "r",
              ##            "r", "r", "r", "r", "r",
              ##            "rw", "rw", "r", "r", "rw", "rw"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)

  if(u$ierr>1) stop("Internal error in 'update.spam.chol.NgPeyton' code ", u$ierr,call.=FALSE)

  if(u$ierr == 1) {
    if (getOption("spam.cholupdatesingular") == "null")
      return(NULL)
    else if (getOption("spam.cholupdatesingular") == "error")
      stop("Singularity problem when updating a Cholesky Factor.")
    else if (getOption("spam.cholupdatesingular") == "warning")
      warning("Singularity problem when updating a Cholesky Factor.\n'object' not updated.")
    else
      stop("'cholupdatesingular' should be 'error', 'null' or 'warning'.")
  }  else {
    slot(object, "entries", check = FALSE) <- u$entries
  }
  invisible(object)
}




solve.spam <- function (a, b,  Rstruct = NULL, ...) {
  nrow <- a@dimension[1]
  ncol <- a@dimension[2]
  if (ncol != nrow)      stop("only square matrices can be inverted")

  if (missing(b)) {
    b <- diag(1, ncol)
  }  else {
    if(!is.matrix(b)) b <- as.matrix(b)
  }
  p <- dim(b)[2]
  if(nrow!=dim(b)[1])stop("'b' must be compatible with 'a'")

  # if we have a spam matrix, we calculate the Cholesky factor
  if (is(a,"spam"))
    if (is(Rstruct, "spam.chol.NgPeyton"))
        a <- update.spam.chol.NgPeyton(Rstruct, a, ...)
    else a <- chol.spam(a, ...)


  if (is(a,"spam.chol.NgPeyton")) {
      # The following is a fast way to perform:
      #     z <- backsolve(a,forwardsolve( t(a),b))
    nsuper <- as.integer( length(a@supernodes)-1)
      if( getOption("spam.force64") || .format.spam(a)$package != "spam" )
          SS <- .format64()
      else
          SS <- .format32

      z <- .C64("backsolves",
                SIGNATURE = c(rep(SS$signature,5),"double", rep(SS$signature, 4),
                              rep("double",3)),

                m = nrow,
                nsuper,
                p,
                a@colindices,
                a@colpointers,
                a@entries,
                a@rowpointers,
                a@invpivot,
                a@pivot,
                a@supernodes,

                vector_dc("double",nrow),
                sol = vector_dc("double",nrow*p),
                as.vector(b,"double"),

                INTENT = c( rep( "r", 10), rep( "rw", 3)),
                NAOK = getOption("spam.NAOK"),
                PACKAGE = SS$package )$sol
  } else z <- backsolve(a, forwardsolve( t(a),b))
    # see the helpfile for a comment about the 't(a)' construct.

  if ( p!=1)    dim(z) <- c(nrow,p)
  return( z)
}

chol2inv.spam <- function (x, ...) {
    ## print("chol2inv.spam")
    nrow <- x@dimension[1]

    if (is(x,"spam.chol.NgPeyton")) {
        y <- vector("double",nrow*nrow)
        y[1L + 0L:(nrow - 1L) * (nrow + 1L)] <- 1.0

        if( getOption("spam.force64") || .format.spam(x)$package !="spam" )
            SS <- .format64()
        else
            SS <- .format32


        z <- .C64("backsolves",
          ##         subroutine backsolves(m,nsuper,nrhs,lindx,xlindx,lnz,
     ## &                   xlnz,invp,perm,xsuper,newrhs,sol,b)
                  SIGNATURE = c( rep( SS$signature, 5), "double", rep( SS$signature, 4),
                                rep( "double", 3)),

                  m = nrow, #r
                  length(x@supernodes)-1, #r
                  nrow, #r
                  x@colindices,
                  x@colpointers,
                  x@entries,
                  x@rowpointers,
                  x@invpivot, #r
                  x@pivot, #r
                  x@supernodes, #r

                  vector_dc("double",nrow), #rw
                  sol = vector_dc("double",nrow*nrow), #w
                  y, #r

                  INTENT = c( rep( "r", 10),
                             "rw", "w", "r"),
                  NAOK = getOption("spam.NAOK"),
                  PACKAGE = SS$package)$sol

        dim(z) <- c(nrow,nrow)
    } else z <- backsolve.spam(x, forwardsolve.spam( t(x), diag(nrow)))
    return( z)
}

backsolve.spam <- function(r, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
# r: spam.chol.NgPeyton structure as returned by chol.spam or a spam object
# x: rhs a vector or a matrix in dense form
# dimensions:  ( m x n) ( n x p)

    if( getOption("spam.force64") || .format.spam(r)$package !="spam" )
        SS <- .format64()
    else
        SS <- .format32


  m <- r@dimension[1]
  if(is.vector(x)) {
    n <- length(x)
    p <- 1
  } else {
    x <- as.matrix(x)
    n <- dim(x)[1]
    p <- dim(x)[2]
  }

  # we separate between "spam.chol.NgPeyton" and "spam"
  if (is(r,"spam.chol.NgPeyton")) {
    if (n!=m) stop("Cholesky factor 'r' not compatible with 'x'")
    nsuper <- length(r@supernodes)-1
    if (!getOption("spam.dopivoting")) {

      z <- .C64("backsolvef",
                SIGNATURE=c(SS$signature, SS$signature, SS$signature, SS$signature,
                            SS$signature, "double", SS$signature, SS$signature,
                            "double"),

                m,
                nsuper,
                p,
                r@colindices,

                r@colpointers,
                r@entries,
                r@rowpointers,
                r@supernodes,

                sol = vector("double",m*p),

                INTENT=c("r", "r", "r", "r", "r", "r", "r", "r", "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE=SS$package)$sol

    }else{
      z <- .C64("pivotbacksolve",
                SIGNATURE=c(SS$signature, SS$signature, SS$signature, SS$signature,
                            SS$signature, "double", SS$signature, SS$signature,
                            SS$signature, SS$signature, "double", "double",
                            "double"),

                m,
                nsuper,
                p,
                r@colindices,

                r@colpointers,
                r@entries,
                r@rowpointers,
                r@invpivot,

                r@pivot,
                r@supernodes,
                vector("double",m),
                sol = vector("double",m*p),

                x,

                INTENT=c("r", "r", "r", "r",
                    "r", "r", "r", "r",
                    "r", "r", "r", "w",
                    "r"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE=SS$package)$sol
    }
  } else {
    if (n!=m) stop("Triangular matrix 'r' not compatible with 'x'")
    # solve R sol = x
    z <- .C64("spamback",
              SIGNATURE=c(SS$signature, SS$signature, "double", "double",
                          "double", SS$signature, SS$signature),
              m=m,
              unused=p,
              sol = vector("double",m*p),
              x=x,

              al=r@entries,
              jal=r@colindices,
              ial=r@rowpointers,

              INTENT=c("rw", "rw", "rw", "rw", "rw", "rw", "rw"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE=SS$package)
    if (z$m<0) stop(gettextf("singular matrix in 'backsolve'. Last zero in diagonal [%d]",
            -z$m), domain = NA)
     else z <- z$sol

  }

  if (p>1) dim(z) <- c(m,p)
  return(z)
}

forwardsolve.spam <- function(l, x,...){#, k = NULL, upper.tri = NULL, transpose = NULL){
#  l: spam.chol.NgPeyton structure as returned by chol.spam
#         or an ordinary lower triangular spam matrix
#  x: rhs a vector a matrix in dense form
#  dimensions:  ( m x n) ( n x p)
#  if (!any(is.null(c(upper.tri,k,transpose ))))
#    warning("'k', 'upper.tri' and 'transpose' argument do not have any effect here")
    if( getOption("spam.force64") || .format.spam(l)$package !="spam" )
        SS <- .format64()
    else
        SS <- .format32

  m <- l@dimension[1]
  if(is.vector(x)) {
    n <- length(x)
    p <- 1L
  } else {
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    p <- dim(x)[2]
  }

  # we separate between "spam.chol.NgPeyton" and "spam"
  if (is(l,"spam.chol.NgPeyton")) {
    if(n!=m) stop("Cholesky factor 'l' not compatible with 'x'")
    nsuper <- length(l@supernodes)-1
    if (!getOption("spam.dopivoting")) {
      z <- .C64("forwardsolvef",
                SIGNATURE=c(SS$signature, SS$signature, SS$signature, SS$signature,
                            SS$signature, "double", SS$signature, SS$signature,
                            "double"),

                m,
                nsuper,
                p,
                l@colindices,

                l@colpointers,
                l@entries,
                l@rowpointers,
                l@supernodes,

                sol = vector("double",m*p),

                INTENT=c("r", "r", "r", "r",
                    "r", "r", "r", "r",
                    "w"),
                NAOK = getOption("spam.NAOK"),
                PACKAGE=SS$package)$sol
    }else{
      z <- .C64("pivotforwardsolve",
                SIGNATURE=c(SS$signature, SS$signature, SS$signature, SS$signature,
                            SS$signature, "double", SS$signature, SS$signature,
                            SS$signature, SS$signature, "double", "double",
                            "double"),
                m,
                nsuper,
                p,
                l@colindices,

                l@colpointers,
                l@entries,
                l@rowpointers,
                l@invpivot,

                l@pivot,
                l@supernodes,
                vector("double",m), #!5#
                sol = vector("double",m*p),

                x,

                INTENT=c("r", "r", "r", "r",
                    "r", "r", "r", "r",
                    "r", "r", "r", "w",
                    "r"),
                PACKAGE=SS$package)$sol
    }
  } else {
    if (n!=m) stop("Triangular matrix 'l' not compatible with 'x'")
    # solve L sol = x
    z <- .C64("spamforward",
              SIGNATURE=c(SS$signature, SS$signature, "double", "double",
                          "double", SS$signature, SS$signature),

              m=m,
              p,
              sol = vector("double",m*p),
              x=x,

              al=l@entries,
              jal=l@colindices,
              ial=l@rowpointers,

              INTENT=c("rw", "r", "w", "r", "r", "r", "r"),
              NAOK = getOption("spam.NAOK"),
              PACKAGE=SS$package)
    if (z$m<0) stop(gettextf("singular matrix in 'forwardsolve'. First zero in diagonal [%d]", -z$m), domain = NA)
    else z <- z$sol
  }
  if (p>1)
    dim(z) <- c(m,p)
  return(z)
}




########################################################################

