\name{kulczynski}
\alias{kulczynski}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Kulczynski distance matrix}
\description{
  Computes Kulczynski distances between the columns of a 0-1-matrix.
}
\usage{
kulczynski(regmat)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{0-1-matrix. Columns are species, rows are regions.}
}
\details{
  The Kulczynski distance between two species is 1-(mean of (number of
  regions where
  both species are present)/(number of regions where species 1 is
  present) and (number of
  regions where
  both species are present)/(number of regions where species 2 is
  present)).
}
\value{
  A symmetrical matrix of Kulczynski distances.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}


\seealso{
  \code{\link{jaccard}}, \code{\link{geco}},\code{\link{qkulczynski}} 
}

\examples{
data(kykladspecreg)
kulczynski(t(kykladspecreg))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line

