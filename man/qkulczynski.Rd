\name{qkulczynski}
\alias{qkulczynski}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Quantitative Kulczynski distance matrix}
\description{
  Computes quantitative Kulczynski distances between the columns of an
  abundance matrix.
}
\usage{
qkulczynski(regmat, log.distance=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{(non-negative) abundance matrix. Columns are species,
    rows are regions.}
  \item{log.distance}{logical. If \code{TRUE}, 1 is added to the
    abundance matrix and then the logs of the values are taken in order
    to compute the distance.}
}
\details{
  The quantitative Kulczynski distance between two species
  is 1-(mean of (mean of over regions minimum abundance of both
  species)/(sum of abundances of species 1) and (mean of over regions
  minimum abundance of both species)/(sum of abundances of species 2)).
  If the abundance matrix is a 0-1-matrix, this gives the standard
  Kulczynski distance.
}
\value{
  A symmetrical matrix of quantitative Kulczynski distances.
}

\references{
  D. P. Faith, P. R. Minchin and L. Belbin (1987) Compositional
  dissimilarity as a robust measure of ecological distance.
  \emph{Vegetation} 69, 57-68.
  }
  
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}


\seealso{
  \code{\link{kulczynski}}
}

\examples{
options(digits=4)
data(kykladspecreg)
qkulczynski(t(kykladspecreg))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line

