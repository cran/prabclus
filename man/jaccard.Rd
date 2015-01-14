\name{jaccard}
\alias{jaccard}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Jaccard distance matrix}
\description{
  Computes Jaccard distances between the columns of a 0-1-matrix.
}
\usage{
jaccard(regmat)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{0-1-matrix. Columns are species, rows are regions.}
}
\details{
  The Jaccard distance between two species is 1-(number of regions where
  both species are present)/(number of regions where at least one
  species is present). As a similarity coefficient, this is S22 in
  Shi (1993).

  Thank you to Laurent Buffat for improving this function!
}
\value{
  A symmetrical matrix of Jaccard distances.
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
  \code{\link{kulczynski}}, \code{\link{dicedist}}
}

\examples{
options(digits=4)
data(kykladspecreg)
jaccard(t(kykladspecreg))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line
