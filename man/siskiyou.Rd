\name{siskiyou}
\alias{siskiyou}
\alias{siskiyou.nb}
\alias{siskiyou.groups}
% \non_function{}
\title{Herbs of the Siskiyou Mountains}
\description{
  Distributions of species of herbs in relation to elevation on quartz
  diorite in the central Siskiyou Mountains. All values are per mille
  frequencies
  in transects (The number of 1 m2 quadrats, among 1000 such quadrats,
  in which a species was observed, based on 1250 1m2 quadrats in the
  first 5  transects, and 400 1m2 quadrats in 6. transect). Observed
  presences in the transect, outside the sampling plots, were coded as
  0.2. Rows correspond to species, columns correspond to regions.
}
\usage{data(siskiyou)}
\format{
  Three objects are generated:
  \item{siskiyou}{numeric matrix giving the 144*6 abundance values.}
  \item{siskiyou.nb}{neighborhood list for the 6 regions.}
  \item{siskiyou.groups}{integer vector of length 144, giving group
    memberships for the 144 species.}
}
\source{
Whittaker, R. H. 1960. Vegetation of the Siskiyou Mountains, Oregon and
California. \emph{Ecol. Monogr}. 30: 279-338 (table 14).
}
\details{
  Reads from example data files \code{LeiMik1.dat, LeiMik1NB.dat,
    LeiMik1G.dat}.
  }

\examples{
data(siskiyou)
}
\keyword{datasets}
