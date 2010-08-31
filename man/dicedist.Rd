\name{dicedist}
\alias{dicedist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Dice distance matrix}
\description{
  Computes a distance derived from Dice's coincidence index
  between the columns of a 0-1-matrix.
}
\usage{
  dicedist(regmat)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{0-1-matrix. Columns are species, rows are regions.}
}
\details{
  The Dice distance between two species is 1 minus the Coincidence
  Index, which is (2*number of regions where
  both species are present)/(2*number of regions where
  both species are present plus number of regions where at least one
  species is present). This is S23 in Shi (1993).
}
\value{
  A symmetrical matrix of Dice distances.
}

\references{
  Shi, G. R. (1993) Multivariate data analysis in palaeoecology and
  palaeobiogeography - a review. \emph{Palaeogeography, 
    Palaeoclimatology, Palaeoecology} 105, 199-234.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{
  \code{\link{kulczynski}},\code{\link{jaccard}}
}

\examples{
data(kykladspecreg)
dicedist(t(kykladspecreg))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line
