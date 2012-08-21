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
  Reads from example data file \code{nb.dat}.
}

\source{
  B. Hausdorf and C. Hennig (2005) The influence of recent geography,
  palaeography and climate on the composition of the faune of the
  central Aegean Islands. \emph{Biological Journal of the Linnean
    Society} 84, 785-795. 
}
\examples{
data(nb)
# nb <- list()
# for (i in 1:34)
#   nb <- c(nb,list(scan(file="(path/)nb.dat",
#                   skip=i-1,nlines=1)))
}
\keyword{datasets}

