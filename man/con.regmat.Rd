\name{con.regmat}
\alias{con.regmat}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Connected regions per species}
\description{
  Returns a vector of the numbers of connected regions per species for a
  presence-absence matrix.   
}
\usage{
con.regmat(regmat, neighbors, count = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{0-1-matrix. Columns are species, rows are regions.}
  \item{neighbors}{A list with a component for every region. The
    components are vectors of integers indicating
    neighboring regions. A region without neighbors (e.g., an island)
    should be assigned a list \code{numeric(0)}.}
  \item{count}{logical. If \code{TRUE}, the number of the currently
    processed species is printed.}
}
\details{
  Uses \code{con.comp}.
}
\value{
  Vector of numbers of connected regions per species.
}

\author{Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/}}

\note{Designed for use in \code{prabtest}.}

\seealso{\code{\link{con.comp}}, \code{\link{prabtest}}}

\examples{
data(nb)
# Note: If you do not use the installed package, replace this by
# nb <- list()
# for (i in 1:34)
#   nb <- c(nb,list(scan(file="(path/)nb.dat",
#                   skip=i-1,nlines=1)))
set.seed(888) 
cp <- cluspop.nb(nb, p.nb=0.1, n.species=10, clus.specs=9,
                 reg.group=1:17,vector.species=c(10))
con.regmat(cp,nb)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
