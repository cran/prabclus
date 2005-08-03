\name{geo2neighbor}
\alias{geo2neighbor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Neighborhood list from geographical distance}
\description{
  Generates a neighborhood list as required by \code{prabinit} from a
  matrix of geographical distances.
}
\usage{
geo2neighbor(geodist,cut=0.1*max(geodist))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{geodist}{\code{dist}-object or symmetric non-negative mmatrix.
    Geographical distances between regions.}
  \item{cut}{non-negative numerical. All pairs of regions with
    \code{distance<=cut} are treated as neighbors.}
}
\value{
  A list of integer vectors, giving the set of neighbors for every region.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\examples{
data(waterdist)
geo2neighbor(waterdist)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line

