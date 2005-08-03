\name{lcomponent}
\alias{lcomponent}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Largest connectivity component}
\description{
  Computes the size of the largest connectivity component of the graph
  of \code{ncol(distmat)} vertices with edges defined by the smallest
  \code{ne} distances.
}
\usage{
lcomponent(distmat, ne = floor(3*ncol(distmat)/4))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{distmat}{symmetric distance matrix.}
  \item{ne}{integer.}
}

\value{
  list with components
  \item{lc}{size of the largest connectivity component.}
  \item{ne}{see above.}
}
\references{
Hennig, C. and Hausdorf, B. (2004) Distance-based parametric bootstrap
tests for clustering of species ranges. \emph{Computational Statistics
  and
  Data Analysis} 45, 875-896.}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{prabtest}}}

\examples{
data(kykladspecreg)
# Note: If you do not use the installed package, replace this by
# kykladspecreg <- read.table("(path/)kykladspecreg.dat")
j <- jaccard(t(kykladspecreg))
lcomponent(j)
}
\keyword{cluster}% at least one, from doc/KEYWORDS




