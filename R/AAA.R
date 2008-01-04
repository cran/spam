# with much input from Roger Bivand for 0.13-2

".onLoad" <-
function (lib, pkg)
{
#  library.dynam("spam",pkg, lib)
	require(methods)
}

".Spam" <- list(eps=.Machine$double.eps,   # smaller than this is considered as zero
                        drop=FALSE,                # drop passed to subset functions
                        printsize=100,             # the max size which we print as regular matrices
                        imagesize=10000,           # the max size which we display as regular matrices
                        trivalues=FALSE,           # with upper./lower/.tri return values (TRUE) or only structure?
                        cex=1200,                  # scaling factor for scatter displays
                        
                        version=list(release=utils:::packageDescription("spam",
                            field="Version"),
			  date=utils:::packageDescription("spam",
                            field="Date")),
                        
                        safemode=TRUE,             # verify double and integer formats. 
                        bcksl=TRUE                 # what type of back/forwardsolve?
                        
                        )

".onAttach" <-
function (lib, pkg)
{
  cat("Package 'spam' is loaded.  Version",
      .Spam$version$release," (", .Spam$version$date, ").",
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

