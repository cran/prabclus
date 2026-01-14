\name{geco}
\alias{geco}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{geco distance matrix}
\description{
  Computes geco distances between the columns of a 0-1-matrix, based on
  a distance matrix between regions (usually, but not necessarily, this is
  a geographical distance).
}
\usage{
geco(regmat,geodist=as.dist(matrix(as.integer(!diag(nrow(regmat))))),
                   transform="piece",
                   tf=0.1,
                   countmode=ncol(regmat)+1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{regmat}{0-1-matrix. Columns are species, rows are regions.}
  \item{geodist}{\code{dist}-object or symmetric non-negative matrix.
    Geographical distances between regions.}
  \item{transform}{transformation applied to the distances before
    computation of geco coefficient, see details. "piece" means
    piecewise linear, namely distance/(\code{tf}*maximum distance)
    if distance<\code{tf}*maximum distance, and 1 otherwise, "log" means
    \code{log((tf*distance)+1)}, "sqrt" means \code{sqrt(tf*distance)},
    "none" means no transformation.}
  \item{tf}{tuning constant for transformation. See \code{transform}.}
  \item{countmode}{optional positive integer. Every 'countmode' algorithm runs
          'geco' shows a message.}
}
\details{
  The geco distance between two species is 0.5*(mean distance
  between region where species 1 is present and closest region where
  species 2 is present plus mean distance
  between region where species 2 is present and closest region where
  species 1 is present). 'closest' to a region could be the regions
  itself. 
  It is recommended (Hennig and Hausdorf, 2006) to transform the
  distances first, because the differences between large distances are
  usually not meaningful or at least much less meaningful than
  differences between small distances for dissimilarity measurement
  between species ranges. See parameter \code{transform}.

  If the between-regions distance is 1 for all pairs of
  non-equal regions, the geco distance degenerates
  to the Kulczynski distance, see \code{kulczynski}.
}
\value{
  A symmetrical matrix of geco distances.
}

\references{
  Hennig, C. and Hausdorf, B. (2006)  A robust distance coefficient
  between distribution areas incorporating geographic distances.
  \emph{Systematic Biology} 55, 170-175.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{
  \code{\link{kulczynski}}
}

\examples{
options(digits=4)
data(kykladspecreg)
data(waterdist)
geco(t(kykladspecreg),waterdist)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line

