\name{prabclus-package}
\alias{prabclus-package}
\alias{prabclus}
\docType{package}
\title{
Test for clustering of presence-absence data
~~ package title ~~
}
\description{
Distance-based parametric bootstrap tests for clustering,
mainly intended for presence-absence data (clustering of
species distribution maps). Jaccard, Kulczynski and geco
distance measures, clustering of presence-absence and abundance
data, and nearest neighbor based noise detection (R port of
Byers and Raftery's (1998) "NNclean"). Main functions are
prabtest (for testing), prabclust (for clustering), prabinit
(for preparing the data) and NNclean (for noise detection). The
help-pages for prabtest and prabclust contain simple standard
executions.
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab prabclus\cr
Version: \tab 2.0-2\cr
Date: \tab 2006-10-06\cr
Depends: \tab R (>= 2.0), MASS\cr
Suggests: \tab mclust02\cr
License: \tab GPL\cr
URL: \tab http://www.homepages.ucl.ac.uk/~ucakche\cr
Packaged: \tab Fri Oct 6 20:57:45 2006; chrish\cr
Built: \tab R 2.6.0; ; 2007-11-05 12:53:48; unix\cr
}

Index:
\preformatted{
NNclean                 Nearest neighbor based clutter/noise detection
autoconst               Spatial autocorrelation parameter estimation
cluspop.nb              Simulation of presence-absence matrices
                        (clustered)
comp.test               Compare species clustering and species groups
con.comp                Connectivity components of an undirected graph
con.regmat              Connected regions per species
distratio               Distance ratio test statistics for distance
                        based clustering
geco                    geco distance matrix
geo2neighbor            Neighborhood list from geographical distance
homogen.test            Classical distance-based test for homogeneity
                        against clustering
hprabclust              Clustering of species ranges from
                        presence-absence matrices (hierarchical
                        methods)
incmatrix               Nestedness matrix
jaccard                 Jaccard distance matrix
kulczynski              Kulczynski distance matrix
kykladspecreg           Snail presence-absence data from Aegean sea
lcomponent              Largest connectivity component
nb                      Neighborhood list for Aegean islands
nbtest                  Test of neighborhood list
nn                      Mean distance to kth nearest neighbor
piecewiselin            Piecewise linear transformation for distance
                        matrices
pop.sim                 p-value simulation for presence-absence
                        matrices clustering test
prabclust               Clustering of species ranges from
                        presence-absence matrices (mixture method)
prabinit                Presence-absence/abundance matrix
                        initialization
prabtest                Parametric bootstrap test for clustering in
                        presence-absence matrices
qkulczynski             Quantitative Kulczynski distance matrix
randpop.nb              Simulation of presence-absence matrices
                        (non-clustered)
summary.prabtest        Print and summary method for prabtest
waterdist               Overwater distances between islands in the
                        Aegean sea
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Christian Hennig <chrish@stats.ucl.ac.uk>

Maintainer: Christian Hennig <chrish@stats.ucl.ac.uk>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
