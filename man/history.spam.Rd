\name{history}
\alias{history.spam}
\alias{spam.history}
\alias{history}
\title{Brief overview of the history}
\description{Brief overview of \code{spam}s  history}
\usage{
  spam.history()
}
\details{This list is maintained for releases 0.05 and higher.
  \item{}{0.15: Induced by need to catch up the devel versions of 0.14.}
  \item{}{0.14: Complete rewrite of \code{chol}, implementation of
    \code{update.spam.chol.NgPeyton} and related routines.
    Rewrite of \code{spam.chol.NgPeyton} class. Proper \code{spam.options}
    handling.  Kronecker product implementation. Rewrite of
    \code{cbind}. Importing and handling of foreign formats.
    Many internal changes for efficieny increases. Major file renaming.
    Some backwards compatibility loss.}
  \item{}{0.13: Minor changes to avoid some warning/errors when
    compiling on different platforms.}
  \item{}{0.12: Many internal changes. Start maintaining a proper
    CHANGELOG file.}
  \item{}{0.11: Considerable help file improvements. Changes in
    \code{c/rbind}, \code{chol}, print/summary methods and other minor
    improvements, first push to CRAN.}
  \item{}{0.10: Minor bug fixes and help file improvements.}
  \item{}{0.09: Proper NAMESPACE handling.}
  \item{}{0.08: \code{fields} uses \code{spam} now as a required
    package.
    The required linking functions are implemented in fields.}
  \item{}{0.07: Change of slot names.}
  \item{}{0.06: Implements \code{spam} and the function linking
            \code{fields} version 3.6. Last stable version before the
	    "merging" with fields. 
              Update to R2.4 and documentation fill-in.}
  \item{}{0.05: Implements \code{spam} and the function linking
              \code{fields} version 3.2. It should represent a fair working
              version, with large gaps in the documentation.}
}
\references{See also a more detailed description at
  \url{http://www.mines.edu/~rfurrer/software/spam/history.shtml}.}
\seealso{\code{\link{todo}}.}
\examples{
\dontrun{
## A complete history of changes
file.show(system.file("NEWS", package = "spam"))
}
}
\keyword{internal}
