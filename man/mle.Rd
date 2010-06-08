% This is file ../spam0.22-0/man/mle.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% written and maintained by Reinhard Furrer.
   

\name{mle}
\alias{neg2loglikelihood.spam}
\alias{mle.spam}
\alias{mle.nomean.spam}
\alias{neg2loglikelihood}
\alias{mle}
\alias{mle.nomean}
\title{Maximum likelihood estimates}
\description{Maximum likelihood estimates of a simple spatial model}
\usage{
neg2loglikelihood.spam(y, X, distmat, Covariance,
                 beta, theta, Rstruct = NULL, ...)

neg2loglikelihood(y, X, distmat, Covariance,
                 beta, theta, ...)

mle.spam(y, X, distmat, Covariance,
     beta0, theta0,
     thetalower, thetaupper, optim.control=NULL,
     Rstruct = NULL,...)

mle(y, X, distmat, Covariance,
     beta0, theta0,
     thetalower, thetaupper, optim.control=NULL, ...)

mle.nomean.spam(y, distmat, Covariance,
     theta0,
     thetalower, thetaupper, optim.control=NULL,
     Rstruct = NULL,...) 

mle.nomean(y, distmat, Covariance,
     theta0,
     thetalower, thetaupper, optim.control=NULL, ...) 
}
\arguments{
  \item{y}{data vector of length n.}
  \item{X}{the design matrix of dimension n x p.}
  \item{distmat}{a distance matrix. Usually the result of a call
     to \code{nearest.dist}.}
  \item{Covariance}{function defining the covariance. See example.}
  \item{beta}{parameters of the trend (fixed effects).}
  \item{theta}{parameters of the covariance structure.}
  \item{Rstruct}{the Cholesky structure of the covariance matrix.}
  \item{beta0,theta0}{inital values.}
  \item{thetalower,thetaupper}{lower and upper bounds of the parameter
    \code{theta}.}
  \item{optim.control}{arguments passed to \code{optim}.}
  \item{...}{additional arguments passed to \code{chol}.}
}
\value{The negative-2-loglikelihood or the output from the function
     \code{optim}.
}
\details{
We provide functions to calculate the
negative-2-log-likelihood and maximum likelihood estimates  for the model

y ~ N_n( X beta, Sigma(h;theta) )

in the case of a sparse or ordinary covariance matrices.

In the case of the \code{*.spam} versions, the covariance function has
to return a \code{spam} object. In the other case, the methods are
correctly overloaded and work either way, slightly slower than the
\code{*.spam} counterparts though.


When working on the sphere, the distance matrix has to be transformed by

h ->   R / 2 sin(h/2)

where R is the radius of the sphere.


The covariance function requires that the first argument is the
distance matrix and the second the parameters, out of which the
first is the range parameter.

In the zero-mean case, the \code{neg2loglikelihood} is calculated by setting
the parameters \code{X} or \code{beta} to zero.
}
%\references{}
\seealso{\code{\link{chol}}
}
\examples{
# True parameter values:
truebeta <- c(1,2,.2)    # beta = (intercept, linear in x, linear in y)
truetheta <- c(.5,2,.02) # theta = (range, sill, nugget)


# Auxiliary covariance function, parameterized by
#    theta = (range, sill, nugget)
spherical <- function(distmat, theta, eps = 1e-06) {
  Sigma <- distmat
  d <- Sigma@entries/theta[1]
 
  Sigma@entries <- ifelse(d < eps,
                     theta[3]+ theta[2],
                   ifelse(d < 1, theta[2]*(1 - 1.5*d + 0.5*d^3), 0))
  return( Sigma)
}

# We now define a grid, distance matrix, and a sample:
x <- seq(0,1,l=5)
locs <- expand.grid( x, x)
X <- as.matrix( cbind(1,locs))  # design matrix
Covariance <- 'spherical'      # covariance function

distmat <- nearest.dist( locs, upper=NULL) # distance matrix
Sigma <- spherical( distmat, truetheta)    # true covariance matrix


set.seed(15)
y <- c(rmvnorm.spam(1,X \%*\% truebeta,Sigma)) # construct sample

# Here is the negative 2 log likelihood:
neg2loglikelihood.spam( y, X, distmat, Covariance,
                       truebeta, truetheta)

# We pass now to the mle:
res <- mle.spam(y, X, distmat, Covariance,
         truebeta, truetheta,thetalower=c(0,0,0),thetaupper=c(1,Inf,Inf))

# Similar parameter estimates here, of course:
mle.nomean.spam(y-X\%*\%res$par[1:3], distmat, Covariance,
         truetheta, thetalower=c(0,0,0), thetaupper=c(1,Inf,Inf))
}
\author{Reinhard Furrer }
\keyword{algebra}

