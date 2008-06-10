### in base:
# dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p=2)
#
### in fields
# rdist( x1, x2)
# rdist.earth(x1, x2, miles = TRUE, R = NULL)
# fields.rdist.near( x1, x2, delta, max.points= NULL)
#
### in sp 
# spDistsN1(pts, pt, longlat=FALSE)
#
### in amap     ### nbproc 	integer, Number of subprocess for parallelization
# Dist(x, method = "euclidean", nbproc = 1, diag = FALSE, upper = FALSE)
#
### in argosfilter
# distance(lat1, lat2, lon1, lon2)  # gc between two pts in km
# distanceTrack(lat,lon)            # gc between pts  in km
#
### in proxy
# dist(x, y = NULL, method = NULL, ..., diag = FALSE, upper = FALSE,
#     pairwise = FALSE, by_rows = TRUE, convert_similarities = TRUE,
#     auto_convert_data_frames = TRUE)
#
### in RFOC
#  GreatDist(LON1, LAT1, LON2, LAT2, EARTHRAD= 6371)
#

nearest.dist <- function( x, y=NULL, method = "euclidean",
             eps = .Spam$eps, delta = 1,
             diag = FALSE, upper = NULL,
             p = 2, miles=TRUE, R=NULL)
{
  # see help for exact parameter meaning

  # diag=FALSE (not include diagonal) only makes sense if y=NULL
  # for the fortran routine we have the following coding:
  # diag=0 include diagonal and y=NULL
  # diag=1 include diagonal and y is given (then the entries are also subject to eps/delta
  # diag=-1 do not include diagonal
  if (!is.na(pmatch(method, "euclidian")))     method <- "euclidean"
  METHODS <- c("euclidean", "maximum", "minkowski", "greatcircle")
  method <- pmatch(method, METHODS)  # result is integer

  if (is.na(method))     stop("invalid distance method")
  if (method == -1)      stop("ambiguous distance method")
  if (delta <= eps)      stop("'delta' should be larger than 'eps'.")
  if (method == 4) {
    if (is.null(R))
      p <- ifelse( miles,3963.34,6378.388)
    else
      p <- R
    if (abs(delta)>180)  stop("'delta' should be smaller than 180 degrees.")
  }
  
  if (is.null(upper)) 
    part <- as.integer(0)
  else
    part <- as.integer( ifelse(upper, 1 ,-1))
  if (is.data.frame(x))  x <- as.matrix(x)
  if (is.list(x)) stop("'x' should be an array or matrix")
           # as.matrix( list() ) does not work
  if (!is.matrix(x)) x <- as.matrix(x)
  nd <- dim(x)[2]
  n1 <- dim(x)[1]

  if (!is.null(y)) {
    # we specify x and y:
    if (is.data.frame(y))  y <- as.matrix(y)
    if (is.list(x)) stop("'x' should be an array or matrix")
    if (!is.matrix(y)) y <- as.matrix(y)
    if (nd!=dim(y)[2]) stop("'x' and 'y' should have the same number of columns.")
    n2 <- dim(y)[1]
    diag <- as.integer(0)
    mi <- min(n1,n2)
    ma <- max(n1,n2)
    nnz <- min(max(.Spam$nearestdistnnz[1],
                   ma*.Spam$nearestdistnnz[2]),
               (as.double(mi)*(mi+1)+(ma-mi)^2)/ ifelse( is.null(upper), 1, 2),
               2^31-2)
    # there is an as.double just in case that mi (and n1 below) is > 2^16
  } else {
    # x = y, i.e. proper distance matrix
    if (n1==1)         stop("More than a single point in 'x' is required.")
    y <- x
    n2 <- n1
    diag <- as.integer( ifelse(diag,-1,1))
    nnz  <- min(max(.Spam$nearestdistnnz[1],
                    n1*.Spam$nearestdistnnz[2]),
                (as.double(n1)*(n1-diag))/ ifelse( is.null(upper), 1, 2),
                2^31-2
                )
  }
  repeat {
    d <- .Fortran("closestdist", nd, as.double(x), n1,  as.double(y), n2, 
                  diag, part,
                  as.double(p[1]), method, 
                  as.double( abs( eps[1])), as.double(abs( delta[1])),
                  colindices=vector("integer",nnz),
                  rowpointers=vector("integer",n1+1),
                  entries=vector("double",nnz),
                  nnz=as.integer(nnz),
                  iflag=as.integer(0),DUP=FALSE,NAOK=!TRUE,
                  PACKAGE="spam")
    
    if (d$iflag==0)
      break
    else {
      if (nnz==2^31-2) stop("distance matrix is too dense (more than 2^31 entries).")
      nnz <- nnz*.Spam$nearestdistincreasefactor*n1/(d$iflag-1)
      warning(paste("You ask for a 'dense' spase distance matrix, I require one more iteration.",
                    "\nTo avoid the iteration, increase the 'nnznearestdist' option\n(constructed ",d$iflag,
                    " lines out of ",n1,").\n",sep=""),
              call. = TRUE)
    }
  }


  dmat <- new("spam")
  slot(dmat,"entries",check=FALSE) <-     d$entries[1:d$nnz]
  slot(dmat,"colindices",check=FALSE) <-  d$colindices[1:d$nnz]
  slot(dmat,"rowpointers",check=FALSE) <- d$rowpointers
  slot(dmat,"dimension",check=FALSE) <-   as.integer(c(n1,n2))
  return( dmat)
}
