\name{summary.prabtest}
\alias{summary.prabtest}
\alias{print.summary.prabtest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Print and summary method for prabtest}
\description{
  Useful compression of the results of \code{prabtest}.
}
\usage{
summary.prabtest(object, ...)

\method{print}{summary.prabtest}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{prabtest}.}
  \item{x}{object of class \code{summary.prabtest}.}
  \item{\dots}{no meaning, necessary for print and summary methods.}
}

\value{
  \code{summary.prabtest} produces an object of class
  \code{summary.prabtest}, which is a list with components
  \item{rrange}{range of the simulation results (test statistic values)
    of \code{object}.}
  \item{rmean}{mean of the simulation results (test statistic values)
    of \code{object}.}
  \item{datac, p.value, pd, tuning, teststat, distance, times, pdfnb}{directly
    taken from \code{object}, see \code{prabtest}.}
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{prabtest}}}

\examples{
data(kykladspecreg)
# Note: If you do not use the installed package, replace this by
# kykladspecreg <- read.table("(path/)kykladspecreg.dat")
data(nb)
# Note: If you do not use the installed package, replace this by
# nb <- list()
# for (i in 1:34)
#   nb <- c(nb,list(scan(file="(path/)nb.dat",
#                   skip=i-1,nlines=1)))
set.seed(5678)
kprab <- prabinit(prabmatrix=kykladspecreg, neighborhood=nb)
kpt <- prabtest(kprab, times=5, pd=0.35, teststat="lcomponent")
# The parameters times and
# pd are only given to fasten the example execution.
summary(kpt)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line


