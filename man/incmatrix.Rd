\name{incmatrix}
\alias{incmatrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Nestedness matrix}
\description{
  Computes species*species nestedness matrix and number of nestings
  (inclusions) from regions*species presence-absence matrix.
}
\usage{
incmatrix(regmat)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{0-1-matrix. Columns are species, rows are regions.}
}

\value{
  A list with components
  \item{m}{0-1-matrix. \code{m[i,j]=1} means that the occupied region of
  species j is a subset (not equal) of the region of species i.}
  \item{ninc}{integer. Number of strict inclusions.}
  \item{neq}{integer. Number of region equalities between species (not
    including equality between species i and i).}
}
\references{
Hausdorf, B. and Hennig, C. (2003) Nestedness of nerth-west European
land snail ranges as a consequence of differential immigration from
Pleistocene glacial refuges. \emph{Oecologia} 135, 102-109.

}
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}



\seealso{
  \code{\link{prabtest}}
}


\examples{
data(kykladspecreg)
incmatrix(t(kykladspecreg))$ninc
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{array}

