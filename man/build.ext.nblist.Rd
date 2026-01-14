\name{build.ext.nblist}
\alias{build.ext.nblist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Internal: generates neighborhood list for diploid loci}
\description{
  This is for use in \code{\link{alleleinit}}.
  Given a neighborhood list of individuals, a new neighborhood list is
  generated in which there are two entries for each individual (entry 1
  and 2 refer to individual one, 3 and 4 to individual 2 and so
  on). Neighborhoods are preserved and additionally the two entries
  belonging to the same individual are marked as neighbors.
}

\usage{
  build.ext.nblist(neighbors,n.individuals=length(neighbors))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{neighbors}{list of integer vectors, where each vector contains
    the neighbors of an individual.}
  \item{n.individuals}{integer. Number of individuals.}
}

  
\value{
  list with \code{2*n.inidividuals} vectors of integers as described
  above.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{\code{\link{alleleinit}}}

\examples{
  data(veronica)
  vnb <- coord2dist(coordmatrix=veronica.coord[1:20,], cut=20,
    file.format="decimal2",neighbors=TRUE)
  build.ext.nblist(vnb$nblist)
}
\keyword{cluster}% at least one, from doc/KEYWORDS




