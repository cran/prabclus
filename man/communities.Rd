\name{communities}
\alias{communities}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Construct communities from individuals}
\description{
  Construct communities from individuals using geographical distance and
  hierarchical clustering. Communities are clusters of geographically
  close individuals, formed by \code{\link{hclust}} with specified
  distance cutoff.
}

\usage{
communities(geodist,grouping=NULL,
                        cutoff=1e-5,method="single")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{geodist}{\code{dist}-object or matrix of geographical distances
    between individuals.}
  \item{grouping}{something that can be coerced into a factor. Different
    groups indicated by \code{grouping} cannot be together in the same
    community. (If \code{NULL}, there is no constraint.)}
  \item{cutoff}{numeric; clustering distance cutoff value, passed on as
    parameter \code{h} to \code{cutree}. Note that if this is smaller
    than the smallest nonzero geographical distance, communities will be
    all sets of individuals that have zero geographical distance to each
    other.}
  \item{method}{\code{method}-parameter for \code{\link{hclust}}.}
}

\value{
  Vector of community memberships for the individuals (integer numbers
  from 1 to the number of communities without interruption.
}

    
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{
  \code{\link{communitydist}}
}

\examples{
  data(veronica)
  ver.geo <- coord2dist(coordmatrix=veronica.coord[1:90,],file.format="decimal2")
  species <-c(rep(1,64),rep(2,17),rep(3,9))
  communities(ver.geo,species)
  
}
\keyword{spatial}% __ONLY ONE__ keyword per line
