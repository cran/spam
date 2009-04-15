# This is file ../spam0.15-4/tests/dist.R
# This file is part of the spam package, 
#      http://www.mines.edu/~rfurrer/software/spam/
# written and maintained by Reinhard Furrer.




options( echo=FALSE)
library( spam, warn.conflict=FALSE)  


distmatrix <- function(x1,x2=NULL)
  {
    if (is.null(x2)) 
      return( as.matrix( dist(x1)))
    else
      return( as.matrix( dist(rbind(x1,x2)))[1:dim(x1)[1],1:dim(x2)[1]+dim(x1)[1]])
  }


test.for.zero <- function( vec, tol = 1.0e-6)
{
  test.value <- sum( (vec^2))
  if( test.value < tol ){
          cat("** PASSED test at tolerance ", tol, fill=TRUE)}
  else{ cat( "## FAILED test value = ", test.value, " at tolerance ", tol,
              fill=TRUE)}

}

########
## as an aside, comparing nearest.dist with dist, use diag=FALSE, upper=TRUE

cat("Results of the form '[1] TRUE' are from 'all.equal'\n\n")

spam.options(printsize=6)


n1 <- as.integer( 4)
n2 <- n1
nd <- as.integer(2)
set.seed(14)
x2 <- x1 <- array(runif(n1*nd), c( n1,nd))

if (F){
# testing the structure
  distmatrix(x1)
  nearest.dist( x1, x1, diag=TRUE,upper=NULL)

# and all other possibilities (2[diag]*3[upper])
# with x1,x1 and x1, NULL:
par(mfcol=c(3,4),mai=c(.15,.15,.05,.05))
display(   nearest.dist( x1, x1,  diag=FALSE, upper=NULL))  # default
display(   nearest.dist( x1, x1,  diag=FALSE, upper=FALSE))
display(   nearest.dist( x1, x1,  diag=FALSE, upper=TRUE))
display(   nearest.dist( x1, x1,  diag=TRUE,  upper=NULL) ) # upper 
display(   nearest.dist( x1, x1,  diag=TRUE,  upper=FALSE))
display(   nearest.dist( x1, x1,  diag=TRUE,  upper=TRUE))
display(   nearest.dist( x1,  diag=FALSE, upper=NULL))
display(   nearest.dist( x1,  diag=FALSE, upper=FALSE))
display(   nearest.dist( x1,  diag=FALSE, upper=TRUE))
display(   nearest.dist( x1,  diag=TRUE,  upper=NULL))
display(   nearest.dist( x1,  diag=TRUE,  upper=FALSE))
display(   nearest.dist( x1,  diag=TRUE,  upper=TRUE))
}


#  nearest.dist( x1) and  nearest.dist( x1,x1) should be identical...
all.equal(  nearest.dist( x1, x1,  diag=FALSE, upper=NULL) ,nearest.dist( x1,  diag=FALSE, upper=NULL) )
all.equal(  nearest.dist( x1, x1,  diag=FALSE, upper=FALSE),nearest.dist( x1,  diag=FALSE, upper=FALSE))
all.equal(  nearest.dist( x1, x1,  diag=FALSE, upper=TRUE) ,nearest.dist( x1,  diag=FALSE, upper=TRUE) )
all.equal(  nearest.dist( x1, x1,  diag=TRUE,  upper=NULL) ,nearest.dist( x1,  diag=TRUE,  upper=NULL) )
all.equal(  nearest.dist( x1, x1,  diag=TRUE,  upper=FALSE),nearest.dist( x1,  diag=TRUE,  upper=FALSE))
all.equal(  nearest.dist( x1, x1,  diag=TRUE,  upper=TRUE) ,nearest.dist( x1,  diag=TRUE,  upper=TRUE) )


# testing  Euclidian
eta <- 1
eps <- 1e-8

o1 <- nearest.dist( x1, diag=TRUE, upper=NULL)
o2 <- distmatrix(x1)

test.for.zero(o2[o2< eta]- o1@entries)

o1 <- nearest.dist( x1, diag=FALSE, upper=NULL)  # is default...
test.for.zero(o2[o2>eps & o2< eta]- o1@entries)

o1 <- nearest.dist( x1, diag=FALSE, upper=!FALSE)  # is default...
o3 <- c(dist(x1))
test.for.zero(o1@entries-o3)



x2 <- x1 <- array(runif(n1*nd), c( n1,nd))
o1 <- nearest.dist( x1,x2, diag=TRUE, upper=NULL)
o2 <- distmatrix(x1,x2)

test.for.zero(o2[o2< eta]- o1@entries)

o1 <- nearest.dist( x1, diag=FALSE, upper=NULL) 
test.for.zero(o2[o2>eps & o2< eta]- o1@entries)

o1 <- nearest.dist( x1, diag=FALSE, upper=!FALSE)  
test.for.zero(o2[o2>eps & o2< eta & lower.tri(o2)]- o1@entries)



# Should cause error:
#    nearest.dist(cbind(1,1))
# this is ok:
test.for.zero( nearest.dist(rbind(1,0)) - c(0,1,0,0))
test.for.zero( nearest.dist(cbind(1,1),cbind(1,0)) -1)



# testing with dist only
test.for.zero( c(as.spam( dist(x1)) - nearest.dist(x1,upper=FALSE,diag=FALSE,delta=2)))


# testing some other norms
method <- "max"
p <- 1.0001
o1 <- nearest.dist( x1,method=method,p=p, diag=FALSE, upper=TRUE )
o3 <- c(dist(x1,method=method,p=p))
test.for.zero(o1@entries-o3)



if (F){ # system.time is not always available...

  n1 <- as.integer( 400)
  set.seed(14)
  x1 <- array(runif(n1*nd), c( n1,nd))
  
  system.time( o1 <- nearest.dist( x1,method="max",p=p) )
  system.time( o1 <- nearest.dist( x1,method="min",p=1) )
  system.time( o1 <- nearest.dist( x1,method="min",p=1.5) )
  system.time( o1 <- nearest.dist( x1,method="min",p=2) )
  system.time( o1 <- nearest.dist( x1,method="euc",p=1) )
  system.time( o1 <- dist( x1) )
}

# testing  GC
n1 <- as.integer( 4)
n2 <- as.integer(6)
set.seed(14)
x1 <- array(runif(n1*2,-90,90), c( n1,2))
x2 <- array(runif(n2*2,-90,90), c( n2,2))



if (F){
# structure
delta <-  180

par(mfcol=c(3,4),mai=c(.15,.15,.05,.05))
eps=0.0001
display( nearest.dist( x1,    eps=eps, delta=delta,method="gr",diag=FALSE,  upper=FALSE))
display( nearest.dist( x1,    eps=eps, delta=delta,method="gr",diag=FALSE,  upper=TRUE))
display( nearest.dist( x1,    eps=eps, delta=delta,method="gr",diag=FALSE,  upper=NULL))
display( nearest.dist( x1,    eps=eps, delta=delta,method="gr",diag=TRUE,   upper=FALSE))
display( nearest.dist( x1,    eps=eps, delta=delta,method="gr",diag=TRUE,   upper=TRUE))
display( nearest.dist( x1,    eps=eps, delta=delta,method="gr",diag=TRUE,   upper=NULL))
display( nearest.dist( x1,x1, eps=eps, delta=delta,method="gr",diag=FALSE,  upper=FALSE))
display( nearest.dist( x1,x1, eps=eps, delta=delta,method="gr",diag=FALSE,  upper=TRUE))
display( nearest.dist( x1,x1, eps=eps, delta=delta,method="gr",diag=FALSE,  upper=NULL))
display( nearest.dist( x1,x1, eps=eps, delta=delta,method="gr",diag=TRUE,   upper=FALSE))
display( nearest.dist( x1,x1, eps=eps, delta=delta,method="gr",diag=TRUE,   upper=TRUE))
display( nearest.dist( x1,x1, eps=eps, delta=delta,method="gr",diag=TRUE,   upper=NULL))
}

# 

if (F){
# if fields would be available, the following can be used as well.
delta <-  180
o2 <- rdist.earth(x1)
o1 <- nearest.dist( x1, method="gr",diag=TRUE,upper=NULL,delta=delta)
test.for.zero(o2- o1@entries)

eps <- (spam.getOption('eps')^.5)*2
o2 <- rdist.earth(x1, R=1)
o1 <- nearest.dist( x1,  method="gr",diag=FALSE,upper=NULL,delta=delta,R=1) 
test.for.zero(o2[o2>eps]- o1@entries[o1@entries>eps])

eps <- 30
o2 <- rdist.earth(x1, R=1)
o1 <- nearest.dist( x1,  method="gr",diag=FALSE,upper=NULL,eps=eps,delta=delta,R=1)
test.for.zero(o2[o2>eps*pi/180]- o1@entries)

delta <- 90
o2 <- rdist.earth(x2,x1,R=1)
o1 <- nearest.dist( x1,x2, method="gr",diag=TRUE,upper=NULL,delta=delta,R=1)
test.for.zero(o2[o2<delta*pi/180]- o1@entries)


}


# correct storage mode conversion
cat("storage conversion:\n")
nx <- 4
x <- expand.grid(as.double(1:nx),as.double(1:nx))
test.for.zero(nearest.dist( x,delta=nx*2,diag=FALSE, upper=TRUE)@entries-
              c(dist(x)))
x <- expand.grid(as.integer(1:nx),as.integer(1:nx))
test.for.zero(nearest.dist( x,delta=nx*2,diag=FALSE, upper=TRUE)@entries-
              c(dist(x)))


# again a bit playing with the parameters:
if (F){
rdist.earth(x1, R=1)
nearest.dist( x1, x1,  method="gr",delta=360,R=1, diag=TRUE)
nearest.dist( x1, x1,  method="gr",delta=360,R=1, diag=FALSE)
nearest.dist( x1,      method="gr",delta=360,R=1, diag=FALSE, upper=FALSE)
nearest.dist( x1,      method="gr",delta=360,R=1, diag=FALSE, upper=TRUE)
nearest.dist( x1,      method="gr",delta=360,R=1, diag=TRUE,  upper=FALSE)
nearest.dist( x1,      method="gr",delta=360,R=1, diag=TRUE,  upper=TRUE)
nearest.dist( x1,x1,   method="gr",delta=360,R=1, diag=FALSE, upper=FALSE)
nearest.dist( x1,x1,   method="gr",delta=360,R=1, diag=FALSE, upper=TRUE)
nearest.dist( x1,x1,   method="gr",delta=360,R=1, diag=TRUE,  upper=FALSE)
nearest.dist( x1,x1,   method="gr",delta=360,R=1, diag=TRUE,  upper=TRUE)
}


# When setting
#     spam.options(nearestdistnnz=c(3,3))
# we do not get a warning because there is a min.
#   nnz <- min(max(.Spam$nearestdistnnz[1], n1 * .Spam$nearestdistnnz[2]),
#              (as.double(n1) * (n1 + 1))/ifelse(is.null(upper),1, 2),
#              2^31 - 2)









#################################################################
# addendum:
# test a few as.spam.dist:

test.for.zero( as.spam(dist(0))@entries -0)
all.equal( as.spam(dist(x1)), nearest.dist(x1,delta=1000))

# diag and upper are only processed when displaying
all.equal( as.spam(dist(x1)), as.spam(dist(x1,diag=TRUE,  upper=TRUE)) )

# 'NA/NaN/Inf's are coerced to zero:
cat("Expect warning:\n")
test.for.zero( as.spam(dist(c(0, NA, 1)))@entries -1)


