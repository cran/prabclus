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
  \item{distmat}{symmetric distance matrix (not a \code{dist}-object).}
  \item{ne}{integer.}
}
\value{
  numerical.
}
\references{
Hennig, C. and Hausdorf, B. (2004) Distance-based parametric bootstrap
tests for clustering of species ranges. \emph{Computational Statistics
  and
  Data Analysis} 45, 875-896.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}


\seealso{\code{\link{prabtest}}}

\examples{
data(kykladspecreg)
# Note: If you do not use the installed package, replace this by
# kykladspecreg <- read.table("(path/)kykladspecreg.dat")
j <- jaccard(t(kykladspecreg))
nn(j,4)
}
\keyword{cluster}% at least one, from doc/KEYWORDS

