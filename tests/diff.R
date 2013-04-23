# This is file ../spam0.29-3/tests/diff.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# written and maintained by Reinhard Furrer.
     








options( echo=FALSE)
library( spam, warn.conflict=FALSE)  

n <- 10
x <- array(rnorm(n^2),c(n,n))

norm(diff(x)-diff(as.spam(x)))
norm(diff(x,d=2)-diff(as.spam(x), d=2))
norm(diff(x,d=4)-diff(as.spam(x), d=4))
norm(diff(x,2, d=2)-diff(as.spam(x),2, d=2))

identical(diff(x,4, d=4), diff(as.spam(x),4, d=4))
