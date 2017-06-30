% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file ../spam/man/toeplitz.Rd
% This file is part of the spam package, 
%      http://www.math.uzh.ch/furrer/software/spam/
% by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
%    Daniel Gerber [ctb], Kaspar Moesinger [ctb]            
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   



\name{validate_spam}
\alias{validate_spam}
\title{Validate a spam objects}
\description{Checks if the spam object has the correct structure.}
\usage{
validate_spam(object)
}
\arguments{
  \item{object}{a spam matrix.}
}
\value{Returns \code{TRUE} if \code{object} is a valid spam objects. }

\seealso{\code{\link{spam.creation}}}
\examples{
validate_spam(spam(1, 20))
}
