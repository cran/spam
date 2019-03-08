# HEADER ####################################################
# This is file spam/R/diff.R.                               #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################

     


########################################################################
diff.spam <- 
function (x, lag = 1, differences = 1, ...) 
{
    xlen <-   dim(x)[1L]
    if (length(lag) > 1L || length(differences) > 1L || lag < 
        1L || differences < 1L) 
        stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen) 
        return( numeric(0))

    for (i in 1L:differences){
      x <- x[(1L+lag):xlen,, drop = FALSE] - x[1L:(xlen-lag),, drop = FALSE] 
      xlen <- xlen - lag
    }
    return( x)
}

setMethod("diff","spam",diff.spam)
