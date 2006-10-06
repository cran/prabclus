\name{pop.sim}
\alias{pop.sim}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{p-value simulation for presence-absence matrices clustering test}
\description{
  Parametric bootstrap simulation of the p-value of a test of a
  homogeneity hypothesis against clustering (or significant nestedness).
  Designed for use within
  \code{prabtest}. The null model is defined by \code{randpop.nb}.
}
\usage{
pop.sim(regmat, neighbors, h0c = 1, times = 200, dist = "kulczynski",
teststat = "isovertice", testc = NULL, geodist=NULL, gtf=0.1,
n.species = ncol(regmat),
specperreg = NULL, regperspec = NULL, species.fixed=FALSE, pdfnb=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{0-1-matrix. Columns are species, rows are regions.}
  \item{neighbors}{A list with a component for every region. The
    components are vectors of integers indicating
    neighboring regions. A region without neighbors (e.g., an island)
    should be assigned a list \code{numeric(0)}.}
  \item{h0c}{numerical. Parameter \code{p.nb} for use in \code{randpop.nb}.}
  \item{times}{integer. Number of simulation runs.}
  \item{dist}{"kulczynski", "jaccard" or "geco", see \code{kulczynski},
    \code{geco} and \code{jaccard}.}
  \item{teststat}{"isovertice", "lcomponent", "distratio", "nn" or
    "inclusions". See
    the corresponding functions, \code{homogen.test} for "isovertice",
    \code{incmatrix} for "inclusions").}
  \item{testc}{numerical. Tuning constant for the test statistics.}
  \item{geodist}{matrix of non-negative reals. Geographical distances
    between regions. Only used if \code{dist="geco"}.}
  \item{gtf}{tuning constant for geco-distance if \code{dist="geco"},
    see \code{"geco"}.}
  \item{n.species}{integer. Number of species.}
  \item{specperreg}{vector of integers. Numbers of species per region
    (is calculated from the data by default).}
  \item{regperspec}{vector of integers. Number of regions per species
    (is calculated from the data by default).}
  \item{species.fixed}{logical. If \code{TRUE}, the sizes of the species
    are taken directly from \code{regmat}. Otherwise, they are drawn by
    random from the empirical distribution of the values from
    \code{regmat}.}
  \item{pdfnb}{logical. Probability correction in \code{randpop.nb}.}
}
\value{
  List with components
  \item{results}{vector of teststatistic values for the simulated matrices.}
  \item{p.above}{p-value if large test statistic leads to rejection.}
  \item{p.below}{p-value if small test statistic leads to rejection.}
  \item{datac}{test statistic value for the original data.}
  \item{testc}{see above.}
}
\references{
Hennig, C. and Hausdorf, B. (2004) Distance-based parametric bootstrap
tests for clustering of species ranges. \emph{Computational Statistics
  and
  Data Analysis} 45, 875-896.
\url{http://stat.ethz.ch/Research-Reports/110.html}.

Hausdorf, B. and Hennig, C. (2003)  Biotic Element Analysis in
Biogeography. \emph{Systematic Biology} 52, 717-723.

Hausdorf, B. and Hennig, C. (2003) Nestedness of north-west European
land snail ranges as a consequence of differential immigration from
Pleistocene glacial refuges. \emph{Oecologia} 135, 102-109.
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}


\seealso{
  \code{\link{prabtest}}, \code{\link{randpop.nb}},
  \code{\link{jaccard}}, \code{\link{kulczynski}},
  \code{\link{homogen.test}}, \code{\link{lcomponent}},
  \code{\link{distratio}}, \code{\link{nn}},
  \code{\link{incmatrix}}.
}

\examples{
data(kykladspecreg)
# Note: If you do not use the installed package, replace this by
# kykladspecreg <- read.table("(path/)kykladspecreg.dat")
data(nb)
# Note: If you do not use the installed package, replace this by
# nb <- list()
# for (i in 1:34)
#   nb <- c(nb,list(scan(file="(path/)nb.dat",
#                   skip=i-1,nlines=1)))
set.seed(1234)
pop.sim(t(kykladspecreg), nb, times=5, h0c=0.35, teststat="nn", testc=3)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
