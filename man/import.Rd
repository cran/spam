% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is file spam/man/import.Rd.                          %
% It is part of the R package spam,                         %
%  --> https://CRAN.R-project.org/package=spam              %
%  --> https://CRAN.R-project.org/package=spam64            %
%  --> https://git.math.uzh.ch/reinhard.furrer/spam         %
% by Reinhard Furrer [aut, cre], Florian Gerber [aut],      %
%    Roman Flury [aut], Daniel Gerber [ctb],                %
%    Kaspar Moesinger [ctb]                                 %
% HEADER END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   










\name{import}
\alias{import}
\alias{read.MM}
\alias{read.HB}
\title{Read External Matrix Formats}
\description{
    Read matrices stored in the Harwell-Boeing or MatrixMarket formats.
}
\usage{
read.HB(file)
read.MM(file)
}
\arguments{
 \item{file}{the name of the file to read, as a character scalar.  

   Alternatively, \code{read.HB} and \code{read.MM} accept connection
   objects.}
}
\value{A sparse matrix of class  \code{spam}.}
\details{The names of files storing matrices in the
  Harwell-Boeing format usually end in \code{".rua"} or \code{".rsa"}.
  Those storing matrices in the MatrixMarket format usually end in
  \code{".mtx"}.

  Currently, only  real assembled Harwell-Boeing can be read with
  \code{read.HB}. Reading MatrixMarket formats is more flexible.
  However, as entries of \code{spam} matrices are of mode \code{double},
  integers matrices are coerced to doubles, patterns lead to matrices
  containing ones and complex are coerced to the real part thereof. In
  these aforementioned cases, a warning is issued.\cr

  MatrixMarket also defines an array format, in which case a (possibly)
  dense \code{spam} object is return (retaining only elements which are
  larger than \code{options('spam.eps')}. A warning is issued. }
\note{
  The functions are based on \code{readHB} and \code{readMM} from
  the library \code{Matrix} to build the connection and read the raw
  data.
  At present, \code{read.MM} is more flexible than  \code{readMM}.
}
\references{
  \url{https://math.nist.gov/MatrixMarket/}
  
  \url{https://sparse.tamu.edu/}
}
\examples{

\dontrun{
image(read.MM(gzcon(url(
  "ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcspwr/bcspwr01.mtx.gz"))))
}

\dontrun{
## Datasets supplied within Matrix
str(read.MM(system.file("external/pores_1.mtx",package = "Matrix")))
str(read.HB(system.file("external/utm300.rua", package = "Matrix")))
str(read.MM(system.file("external/lund_a.mtx", package = "Matrix")))
str(read.HB(system.file("external/lund_a.rsa", package = "Matrix")))
}
}
\author{Reinhard Furrer based on \code{Matrix} functions by
  Douglas Bates \email{bates@stat.wisc.edu} and Martin Maechler \email{maechler@stat.math.ethz.ch}}
\keyword{IO}
\keyword{array}
\keyword{algebra}

