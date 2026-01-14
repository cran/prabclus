\name{nbtest}
\alias{nbtest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Test of neighborhood list}
\description{
  Tests a list of neighboring regions for proper format. Neighborhood is
  tested for being symmetrical. Causes an error if tests fail.
}
\usage{
nbtest(nblist, n.regions=length(nblist))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{nblist}{A list with a component for
    every region. The
    components are vectors of integers indicating
    neighboring regions. A region without neighbors (e.g., an island)
    should be assigned a vector \code{numeric(0)}.}
  \item{n.regions}{Number of regions.}
}

\value{
  \code{invisible{TRUE}}.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{\code{\link{prabinit}}.}

\examples{
data(nb)
nbtest(nb)
nb[[1]][1] <- 1
try(nbtest(nb))
}
\keyword{spatial}% at least one, from doc/KEYWORDS
