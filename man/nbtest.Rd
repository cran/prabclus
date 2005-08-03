\name{nbtest}
\alias{nbtest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Test of neighborhood list}
\description{
  Tests a list of neighboring regions for proper format. Neighborhood is
  tested for being symmetrical. Causes an error if tests fail.
}
\usage{
nbtest(nblist, n.regions=length(nblist))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{nblist}{A list with a component for
    every region. The
    components are vectors of integers indicating
    neighboring regions. A region without neighbors (e.g., an island)
    should be assigned a vector \code{numeric(0)}.}
  \item{n.regions}{Number of regions.}
}

\value{
  \code{invisible{TRUE}}.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{prabinit}}.}

\examples{
data(nb)
# Note: If you do not use the installed package, replace this by
# nb <- list()
# for (i in 1:34)
#   nb <- c(nb,list(scan(file="(path/)nb.dat",
#                   skip=i-1,nlines=1)))
nbtest(nb)
nb[[1]][1] <- 1
try(nbtest(nb))
}
\keyword{spatial}% at least one, from doc/KEYWORDS
