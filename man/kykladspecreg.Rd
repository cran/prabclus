\name{kykladspecreg}
\alias{kykladspecreg}
% \non_function{}
\title{Snail presence-absence data from Aegean sea}
\description{
  0-1-matrix where rows are snail species and columns are islands in the
  Aegean sea. An entry of 1 means that the species is present in the
  region.     
}
\usage{data(kykladspecreg)}
\format{
  A 0-1 matrix with 80 rows and 34 columns.}
\source{
  B. Hausdorf and C. Hennig (2005) The influence of recent geography,
  palaeography and climate on the composition of the faune of the
  central Aegean Islands. \emph{Biological Journal of the Linnean
    Society} 84, 785-795.
}
\details{
  Reads from example data file \code{kykladspecreg.dat}.
  }

\seealso{
  \code{\link{nb}} provides neighborhood information about the 34
  islands.
  }
\examples{
data(kykladspecreg)
}
\keyword{datasets}
% \keyword{spatial}

