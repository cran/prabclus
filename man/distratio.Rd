\name{distratio}
\alias{distratio}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Distance ratio test statistics for distance based clustering}
\description{
  Calculates the ratio between the \code{prop} smallest and largest
  distances of a distance matrix.
}
\usage{
distratio(distmat, prop = 0.25)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{distmat}{symmetric distance matrix.}
  \item{prop}{numerical. Proportion between 0 and 1.}
}

\details{
  Rounding is by \code{floor} for small and \code{ceiling} for large
  distances.
}

\value{
  A list with components
  \item{dr}{ratio of \code{prop} smallest to \code{prop} largest distances.}
  \item{lowmean}{mean of \code{prop} smallest distances.}
  \item{himean}{mean of \code{prop} smallest distances.}
  \item{prop}{see above.}
}
\references{
Hennig, C. and Hausdorf, B. (2002) Distance-based parametric bootstrap
tests for clustering of species ranges, submitted,
\url{http://stat.ethz.ch/Research-Reports/110.html}.}

\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\seealso{\code{\link{prabtest}}}

\examples{
data(kykladspecreg)
# Note: If you do not use the installed package, replace this by
# kykladspecreg <- read.table("(path/)kykladspecreg.dat")
j <- jaccard(t(kykladspecreg))
distratio(j)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

