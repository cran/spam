# Framework introduce with much input from Roger Bivand for 0.13-2 and higher.

".onLoad" <- function (lib, pkg) {
        require(methods) 
}


spam.Version <- function() {
  release <- utils:::packageDescription("spam",field="Version")
  date <- utils:::packageDescription("spam",field="Date")
  list(status="",
              major=sub("-","",substr(release,1,4)),
              minor=substr(sub("-","",substr(release,5,7)),1,1),
              year=substr(date,1,4),
              month=substr(sub("200.-","",date),1,2),
              day=sub("200.-..-","",date),
              version.string= paste("Spam version ",
                utils:::packageDescription("spam",field="Version")," (",
                utils:::packageDescription("spam",field="Date"),")",sep='')
              )
}
                
spam.version <- spam.Version()
class(spam.version) <- "simple.list"
  
".Spam" <- list(eps=.Machine$double.eps,   # smaller than this is considered as zero
                drop=FALSE,                # drop passed to subset functions
                printsize=100,             # the max size which we print as regular matrices
                imagesize=10000,           # the max size which we display as regular matrices
                trivalues=FALSE,           # with upper./lower/.tri return values (TRUE) or only structure?
                cex=1200,                  # scaling factor for scatter displays
                
             
                safemode=c(TRUE,TRUE,TRUE),  # verify double and integer formats and else...
                dopivoting=TRUE,           # what type of back/forwardsolve?
                cholsymmetrycheck=TRUE,     # Should symmetry be tested in the cholesky factorization
                cholpivotcheck=TRUE,        # Should the pivot be tested?
                cholupdatesingular="warning",     # ("error", "warning","NULL")
                cholincreasefactor=c(1.25,1.25),
                nearestdistincreasefactor=1.25,
                nearestdistnnz=c(400^2,400)
                )
#noquote(unlist(format(.Spam[-1])) )


".onAttach" <- function (lib, pkg) {
   cat("Package 'spam' is loaded. ", spam.version$version.string,".",
       "\nType demo( spam) for some demos,",
       " help( Spam) for an overview of this package.\n",
       sep='',fill=TRUE)
   unlockBinding(".Spam", asNamespace("spam"))
 }

"spam.getOption" <- function (x)
  spam.options(x)[[1]]


"spam.options" <- function (...) {
    if (nargs() == 0) return(.Spam)
    current <- .Spam
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.Spam[arg]),
               stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
#    changed <- current[n]  #rf
    current[n] <- temp
    if (sys.parent() == 0) env <- asNamespace("spam") else env <- parent.frame()
    assign(".Spam", current, envir = env)
    invisible(current)
}

#      library(spam,lib='~/todelete/0.15R7')
#      .Spam=spam.options()
int0 <- as.integer(0)
int1 <- as.integer(1)
int2 <- as.integer(2)
