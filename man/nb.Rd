\name{nb}
\alias{nb}
% \non_function{}
\title{Neighborhood list for Aegean islands}
\description{
  List of neighboring islands for 34 Aegean islands.
}
\usage{data(nb)}
\format{
  List with 34 components, all being vetors of integers (or
  \code{numeric(0)} in case of no neighbors) indicating the neighboring
  islands.
}
\details{
  Reads from example data file \code{kyklneigh.txt}.
}

\source{
  Hausdorf, B. and Hennig, C. (2002) Biogeographical tests of the
  Importance of Vicariance in the Diversification of Mediterranean Land
  Snails. Submitted. 
  
}
\examples{
data(nb)
# Note: If you do not use the installed package, replace this by
# nb <- list()
# for (i in 1:34)
#   nb <- c(nb,list(scan(file="(path/)nb.dat",
#                   skip=i-1,nlines=1)))
}
\keyword{datasets}

