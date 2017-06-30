# HEADER ####################################################
# This file is part of the spam package,                    #
#      http://www.math.uzh.ch/furrer/software/spam/         #
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      #
#    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            #
# HEADER END ################################################


validspamobject <- function( ...) {
#    .Deprecated('validate_spam()')
    message("`validspamobject()` is deprecated. Use `validate_spam()` directly")
    validate_spam( ...)
    }

spam.getOption <- function(...) {
#    .Deprecated(msg="`spam.getOption( arg)` is deprecated.\n Use `getOption( spam.arg)` directly")
    message("`spam.getOption( arg)` is deprecated. Use `getOption( spam.arg)` directly")
    getOption(...)
    
}
spam.options <- function(...) {
#    .Deprecated(msg="`spam.options( arg)` is deprecated.\n Use `options( spam.arg)` directly")
    message("`spam.options( arg)` is deprecated. Use `options( spam.arg)` directly")

    if (nargs() == 0) {
        tmp <- names(options())
        return(options()[ grep( "spam.", tmp)])
    }

    args <- list(...)
    names(args) <- paste0("spam.", names(args))
    do.call("options", args)
}    

