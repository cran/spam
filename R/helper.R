# This is file ../spam0.15-5/R/helper.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# written and maintained by Reinhard Furrer.







########################################################################
########################################################################
# a few nice helper functions:



adiag.spam <- function(...){
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
    tmp <- adiag.spam( args[[1]], args[[2]])
    for ( i in 3:nargs)
      tmp <- adiag.spam( tmp, args[[i]])
    return( tmp)
  }
}

# Draw from a multivariate normal:
# (Algorithm 2.3 from Rue and Held)
rmvnorm.spam <- function(n,mu=rep(0, nrow(Sigma)), Sigma,...) {
  # taken from ?chol.spam
  cholS <- chol.spam( Sigma,...)
  # cholS is the upper triangular part of the permutated matrix Sigma
  iord <- ordering(cholS, inv=TRUE)

  N <- dim(Sigma)[1]

  R <- as.spam(cholS)
  retval <- ( array(rnorm(n*N),c(n,N)) %*% R)[,iord,drop=F]
     # It is often better to order the sample than the matrix
     # R itself.
  return(sweep(retval, 2, mu, "+"))
}

# Draw from a multivariate normal given a precision matrix:
# (Algorithm 2.3 from Rue and Held)
rmvnorm.prec <- function(n,mu=rep(0, nrow(Q)), Q, ...) {

  Lt <- chol( Q,...)
  # Lt is the upper triangular part of the permutated matrix Sigma

  N <- dim(Q)[1]
  nu <- backsolve(Lt, array(rnorm(n*N),c(N,n)))
  return(t(nu+mu))
}


# Draw from the canonical representation of a multivariate normal:
# (Algorithm 2.5 from Rue and Held)
rmvnorm.canonical <- function(n, b, Q, ...) {
  N=dim(Q)[1]
  Lt <- chol(Q,...)
  if(is(Lt,"spam.chol.NgPeyton")) {
     mu <- drop(solve(Lt,b))	
  } else {
     mu <- backsolve(Lt,forwardsolve( t(Lt),b))
  }    
  nu <- backsolve(Lt, array(rnorm(n*N),c(N,n)))
  return(t(nu+mu))
}
