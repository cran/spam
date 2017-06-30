# HEADER ####################################################
# This is file  spam/R/rep_len_long.R.                      #
# This file is part of the spam package,                    #
#      http://www.math.uzh.ch/furrer/software/spam/         #
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb],      #
#    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            #
# HEADER END ################################################


rep_len_long <- function(x, length.out, NAOK = getOption("spam.NAOK")){
    if(length.out <= 2147483647){
        return(rep_len(x = x, length.out = length.out))
    }
    return(
        .C64("rep_len_long",
             SIGNATURE = c("double", "int64", "int64",
                 "double"),
             
             x = x,
             lx = length(x),
             length = length.out,
             
             out = numeric_dc(length.out),
             
             INTENT = c("r","r","r",
                 "w"),
             NAOK = NAOK,
             PACKAGE = "spam"
             )$out
        )
    
}
