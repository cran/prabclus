\name{nn}
\alias{nn}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mean distance to kth nearest neighbor}
\description{
  Computes the mean of the distances from each point to its \code{ne}th
  nearest neighbor.
}
\usage{
nn(distmat, ne = 1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{distmat}{symmetric distance matrix.}
  \item{ne}{integer.}
}
\value{
  numerical.
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
nn(j,4)
}
\keyword{cluster}% at least one, from doc/KEYWORDS

