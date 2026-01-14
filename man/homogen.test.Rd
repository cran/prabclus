\name{homogen.test}
\alias{homogen.test}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Classical distance-based test for homogeneity against clustering}
\description{
  Classical distance-based test for homogeneity against clustering. Test
  statistics is number of isolated vertices in the graph of smallest
  distances. The homogeneity model is a random graph model where \code{ne}
  edges are drawn from all possible edges.
}
\usage{
homogen.test(distmat, ne = ncol(distmat), testdist = "erdos")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{distmat}{numeric symmetric distance matrix.}
  \item{ne}{integer. Number of edges in the data graph, corresponding to
    smallest distances.}
  \item{testdist}{string. If \code{testdist="erdos"}, the test
    distribution is a Poisson asymptotic distibution as given by Erdos
    and Renyi (1960). If \code{testdist="ling"}, the test distribution
    is exact as given by Ling (1973), which needs much more computing time.}
}

\details{
  The "ling"-test is one-sided (rejection if the number of isolated vertices
  is too large), the "erdos"-test computes a one-sided as well as a
  two-sided p-value.
}

\value{
  A list with components
  \item{p}{p-value for one-sided test.}
  \item{p.twoside}{p-value for two-sided test, only if \code{testdist="erdos"}.}
  \item{iv}{number of isolated vertices in the data.}
  \item{lambda}{parameter of the Poisson test distribution, only if
    \code{testdist="erdos"}.}
  \item{distcut}{largest distance value for which an edge has been
    drawn.}
  \item{ne}{see above.}
}
\references{
  Erdos, P. and Renyi, A. (1960) On the evolution of random
  graphs. \emph{Publications of the Mathematical Institute of the
  Hungarian Academy of Sciences} 5, 17-61.

  Godehardt, E. and Horsch, A. (1995) Graph-Theoretic Models for Testing
  the Homogeneity of Data. In Gaul, W. and Pfeifer, D. (Eds.) \emph{From
  Data to Knowledge}, Springer, Berlin, 167-176.
  
  Ling, R. F. (1973) A probability theory of cluster
  analysis. \emph{Journal of the American Statistical Association} 68,
  159-164.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}


\seealso{\code{\link{prabtest}}}

\examples{
options(digits=4)
data(kykladspecreg)
j <- jaccard(t(kykladspecreg))
homogen.test(j, testdist="erdos")
homogen.test(j, testdist="ling")
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line







