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
  Hausdorf, B. and Hennig, C. (2002) Biogeographical tests of the
  Importance of Vicariance in the Diversification of Mediterranean Land
  Snails. Submitted. 
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
# Note: If you do not use the installed package, replace this by
# kykladspecreg <- read.table("(path/)kykladspecreg.dat")
}
\keyword{datasets}
% \keyword{spatial}

