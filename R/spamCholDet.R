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


# The main function here is very similar to chol.spam().
# We work with different files such that changes can be easily ported via `meld`.

determinant.spam <- function(x, logarithm = TRUE, pivot = "MMD",method="NgPeyton",
                              memory=list(),eps = getOption("spam.eps"), ...){
   # print("determinant.spam")

  if (eps<.Machine$double.eps) stop("'eps' should not be smaller than machine precision",call.=FALSE)
  logdet <- list()
 #### start from above
  nrow <- x@dimension[1]
  nnzA <- as.integer( x@rowpointers[nrow+1]-1)
  if(nrow!=x@dimension[2]) stop("non-square matrix in 'chol'",call.=FALSE)

  if(nrow <= 1) {
    return(determinant(as.matrix(x@entries))) }

  if(getOption("spam.cholsymmetrycheck")) {
    test <- isSymmetric.spam(x, tol = eps*100)
    if (!isTRUE(test))
      stop("Input matrix to 'chol' not symmetric (up to 100*eps)",call.=FALSE)
  }

  if (method != "NgPeyton")
    warning(gettextf("method = '%s' is not supported. Using 'NgPeyton'",
                     method), domain = NA)

  if (length(pivot)==nrow) {
    doperm <- 0L
    pivot <- as.vector(pivot,"integer")
    if (getOption("spam.cholpivotcheck")) {
      checkpivot(pivot,nrow)
    }
  } else if (length(pivot)==1) {
    if (pivot==FALSE) {
      doperm <- 0L
      pivot <- seq_len(nrow)
    } else if(pivot==TRUE) {
      doperm <- 1L
      pivot <- vector("integer",nrow)
    } else {
      doperm <- as.integer( switch(match.arg(pivot,c("MMD","RCM")),MMD=1,RCM=2))
      pivot <- vector("integer",nrow)
    }
  } else stop("'pivot' should be 'MMD', 'RCM' or a permutation")


  #!6#
    nnzcfact <- c(5,1,5)
    nnzRfact <- c(5,1,2)

  # nnzcolindices = length of array holding the colindices
  if(is.null(memory$nnzcolindices))  {
    nnzcolindices <- ifelse((nnzA/nrow < 5), # very sparse matrix
                            max(1000,nnzA*(1.05*nnzA/nrow-3.8)),
                            nnzA)*nnzcfact[doperm+1]
    nnzcolindices <- max(nnzcolindices,nnzA)
 }else {
    nnzcolindices <- max(memory$nnzcolindices,nnzA)
    memory$nnzcolindices <- NULL
  }
  # nnzR = length of array holding the nonzero values of the factor
  if(is.null(memory$nnzR))    nnzR <- min(max(4*nnzA,floor(.4*nnzA^1.2))*nnzRfact[doperm+1],nrow*(nrow+1)/2)  else {
    nnzR <- memory$nnzR
    memory$nnzR <- NULL
  }
  if(is.null(memory$cache))    cache <- 64  else {
    cache <- memory$cache
    memory$cache <- NULL
  }

  if (length( memory)>0 )
    warning("The component(s) ", paste("'",names(memory),"'",sep='',collapse=","),
            " of the argument 'memory'\npassed to function 'chol' not meaningful and hence ignored.",call.=FALSE)
  ## print("determinant.spam")
  ## z <- .Fortran("cholstepwise",
  ##               nrow = as.integer(nrow) ,nnzA = as.integer(x@rowpointers[nrow+1]-1),
  ##               d =  as.double(x@entries),jd = as.integer(x@colindices),id = as.integer(x@rowpointers),
  ##               doperm = as.integer(doperm), invp = vector("integer",nrow), perm = as.integer(pivot),
  ##               nnzlindx = vector("integer",1),
  ##               nnzcolindices = as.integer(nnzcolindices),
  ##               lindx = vector("integer",nnzcolindices),
  ##               xlindx = vector("integer",nrow+1),     #
  ##               nsuper = vector("integer",1),          #
  ##               nnzR = as.integer(nnzR),#
  ##               lnz = vector("double",nnzR),        #
  ##               xlnz = vector("integer",nrow+1),     #
  ##               snode = vector("integer",nrow),
  ##               xsuper = vector("integer",nrow+1),
  ##               cachesize = as.integer(cache),
  ##               ierr = 0L,
    ##               NAOK = getOption("spam.NAOK"), PACKAGE = "spam")

    if( getOption("spam.force64") || prod(dim(x)) > 2147483647)
        SS <- .format64()
    else
        SS <- .format32
    cholstep.intent <- c("r", "r", "r", "r",
                    "r", "rw", "rw", "rw",
                    "rw", "rw", "rw", "rw",
                    "rw", "rw", "rw", "rw",
                    "rw", "rw", "rw", "rw")
    cholstep.signature <- rep(SS$signature, 20)
    cholstep.signature[c(3,15)] <- "double"

    z <- .C64("cholstepwise",
      ##             subroutine cholstepwise(m,nnzd,
     ## &     d,jd,id,    doperm,invp,perm,
     ## &                nsub,nsubmax,
     ## &                lindx,xlindx,nsuper,nnzlmax,lnz,xlnz,
     ## &                snode,xsuper,
     ## &                cachsz,ierr)
              SIGNATURE = cholstep.signature,

              nrow = nrow ,
              nnzA = x@rowpointers[nrow+1]-1,
              d =  x@entries,
              jd = x@colindices,

              id = x@rowpointers,
              doperm = doperm,
              invp = vector_dc( SS$type, nrow),
              perm = pivot,

              nnzlindx = vector_dc( SS$type, 1),
              nnzcolindices = nnzcolindices,
              lindx = vector_dc( SS$type, nnzcolindices),
              xlindx = vector_dc( SS$type,nrow+1),     #

              nsuper = vector_dc( SS$type,1),          #
              nnzR = nnzR,#
              lnz = vector_dc( "double", nnzR),        #
              xlnz = vector_dc( SS$type, nrow+1),     #

              snode = vector_dc( SS$type, nrow),
              xsuper = vector_dc( SS$type, nrow+1),
              cachesize = cache,
              ierr = 0,

              ## INTENT = cholstep.intent,
              NAOK = getOption("spam.NAOK"),
              PACKAGE = SS$package)


  if(z$ierr == 1) stop("Singularity problem when calculating the Cholesky factor.")
  if(z$ierr == 6) stop("Inconsitency in the input",call.=FALSE)

  while( z$ierr>1) {
    if(z$ierr == 4) {
      warning("Increased 'nnzR' with 'NgPeyton' method\n",
              "(currently set to ",nnzR," from ",ceiling(nnzR*getOption("spam.cholpar")[1]),")",call.=FALSE)
      nnzR <- ceiling(nnzR*getOption("spam.nnzRinc"))
    }
    if(z$ierr == 5) {
      warning("Increased 'nnzcolindices' with 'NgPeyton' method\n",
         "(currently set to ",nnzcolindices," from ",ceiling(nnzcolindices*getOption("spam$cholpar")[2]),")",call.=FALSE)
      nnzcolindices <- ceiling(nnzcolindices*getOption("spam.cholpar")[2])
    }
    print("whileloop in determinant.spam") ##TODO find a case with z$ierr > 1 before migration possible
    z <- .Fortran("cholstepwise",
                  nrow = nrow,nnzA = as.integer(x@rowpointers[nrow+1]-1),
                  d =  as.double(x@entries),jd = x@colindices,id = x@rowpointers,
                  doperm = doperm,invp = vector("integer",nrow), perm = pivot,
                  nnzlindx = vector("integer",1),
                  nnzcolindices = as.integer(nnzcolindices),
                  lindx = vector("integer",nnzcolindices),
                  xlindx = vector("integer",nrow+1),     #
                  nsuper = vector("integer",1),          #
                  nnzR = as.integer(nnzR),#
                  lnz = vector("double",nnzR),        #
                  xlnz = vector("integer",nrow+1),     #
                  snode = vector("integer",nrow),
                  xsuper = vector("integer",nrow+1),
                  cachesize = as.integer(cache),
                  ierr = 0L,
                  NAOK = getOption("spam.NAOK"), PACKAGE = "spam")

  }
 #### end from above
  if(z$ierr == 1) {
                                        # all other errors trapped
      warning("singularity problem or matrix not positive definite",call.=FALSE)
      logdet$modulus <- NA
   } else{
    tmp <- 2* sum( log( z$lnz[ z$xlnz[ -(z$nrow+1)]]))
    if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)
  }

  attr(logdet$modulus,"logarithm") <- logarithm

  logdet$sign <- ifelse(z$ierr == 1,NA,1)
  attr(logdet,"class") <- "det"

  return(logdet)
}

determinant.spam.chol.NgPeyton <- function(x, logarithm = TRUE,...)
{
  logdet <- list()


  tmp <- sum( log(x@entries[ x@rowpointers[-(x@dimension[1]+1)]]))
  if (logarithm) logdet$modulus <- tmp else logdet$modulus <- exp(tmp)

  attr(logdet$modulus,"logarithm") <- logarithm

  logdet$sign <- 1
  attr(logdet,"class") <- "det"

  return(logdet)
}

setMethod("determinant","spam",               determinant.spam)
setMethod("determinant","spam.chol.NgPeyton", determinant.spam.chol.NgPeyton)

## The ``Right Thing'' to do :
## base::det() calls [base::]determinant();
## our det() should call our determinant() :
det <- base::det
environment(det) <- environment()## == asNamespace("Matrix")
######################################################################
########################################################################

