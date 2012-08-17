# This is file ../spam0.29-2/R/makeprec.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# written and maintained by Reinhard Furrer.
     



"make.prec" <- function(par, dims, model = "m1p1",  eps = .Spam$eps){
  if (model=="m1p1"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    y <- numeric(prod(dims))
    y[dims[1]+1] <- -par[1] 
    return( kronecker(diag.spam(dims[2]), toeplitz.spam(x,eps=eps))+toeplitz.spam(y,eps=eps))
  }
  
  if (model=="m1p2"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    y <- numeric(prod(dims))
    y[dims[1]+1] <- -par[2] 
    return( kronecker(diag.spam(dims[2]), toeplitz.spam(x,eps=eps))+toeplitz.spam(y,eps=eps))
  }
  
  if (model=="m2p3"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    
    y <- numeric(dims[1])
    y[1:2] <- c(-par[2],-par[3])

    z <- numeric(dims[2])
    z[2] <- 1


    p1 <- kronecker( diag.spam(dims[2]), toeplitz.spam(x,eps=eps))
    p2 <- kronecker( toeplitz.spam(z,eps=eps), toeplitz.spam(y,eps=eps))
    return( p1 + p2)
  }
  if (model=="m2p4"){
    x <- numeric(dims[1])
    x[1:2] <- c(1,-par[1])
    
    y <- numeric(dims[1])
    y[1:2] <- c(-par[2],-par[3])
    w <- numeric(dims[1])
    w[1:2] <- c(-par[2],-par[4])

    z <- numeric(dims[2])
    z[2] <- 1


    p1 <- kronecker( diag.spam(dims[2]), toeplitz.spam(x,eps=eps))
    p2 <- kronecker( toeplitz.spam(z,rep(0,dims[2]),eps=eps), toeplitz.spam(y,w,eps=eps))
    p3 <- kronecker( toeplitz.spam(rep(0,dims[2]),z,eps=eps), toeplitz.spam(w,y,eps=eps))
    return( p1 + p2 + p3)
  }
  stop("Model not implemented yet!")

}
