\name{phipt}
\alias{phipt}
\alias{cfchord}
\alias{shared.problist}
\alias{diploidcomlist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Distances between communities, auxiliary functions}
\description{
  Auxiliary functions for \code{\link{communitydist}}. \code{phipt}
  computes phiPT/phiST (Peakall and Smouse, 2012, Meirmans,
  2006) between two communities. \code{cfchord} computes the
  chord-distance (Cavalli-Sforza and Edwards, 1967) between two lists or
  locus-wise relative allele frequencies. \code{shared.problist}
  computes a straightforward generalisation of the shared allele
  distance (Bowcock et al., 1994) between
  individuals for communities, namely the `overlap', i.e., sum of the
  minima of the
  allele relative frequencies. \code{diploidcomlist} constructs the
  input lists for \code{cfchord} and \code{shared.problist} from an
  \code{alleleobject}. It provides relative frequencies for all
  alles of all loci in all communities. 
}

\usage{

phipt(alleleobj,comvector,i,j)
cfchord(p1,p2)
shared.problist(p1,p2)
diploidcomlist(alleleobj,comvector,diploid=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{alleleobj}{if \code{diploid=TRUE}, an object of class
    \code{alleleobject} as produced by
    function  \code{\link[prabclus]{alleleinit}}. This has the required
    information on the individuals that are grouped into communities. In
    case \code{diploid=FALSE}, a list that needs to have components
    \code{n.variables} (number of loci), \code{alevels} (vector of
    allele names, see \code{\link{alleleinit}}) and
    \code{charmatrix} (matrix of characters with one row for every
    individual and one column for every locus giving the alleles; see
    examples below for how this can be constructed for a
    \code{prabobject} with presence-absence data).}
  \item{comvector}{vector of integers indicating to which
    community an individual belongs.}
  \item{i}{integer. Number of community.}
  \item{j}{integer. Number of community. The phiPT-distance is computed
    between the communities numbered \code{i} and \code{j}}
  \item{p1}{list. Every list entry refers to a locus and is a vector of
    relative frequencies of the alleles present in that locus in a
    community.} 
  \item{p2}{list. Every list entry refers to a locus and is a vector of
    relative frequencies of the alleles present in that locus in a
    community. The chord or shared allele distance is computed between
    the communities encoded by \code{p1} and \code{p2}.} 
  \item{diploid}{logical, indicating whether loci are diploid, see
    \code{alleleobj}.}

}

  

\value{
  \code{cfchord} gives out the value of the chord
  distance. \code{shared.problist} gives out the distance
  value. \code{diploidcomlist} gives out a two-dimensional list. The
  list has one entry for each community, which is itself a list. This
  community list has one entry for each locus, which is a vector that
  gives the relative frequencies of the different alleles in    
  \code{phipt} gives out a list with components \code{phipt, vap, n0,
  sst, ssg, msa, msw}. These refer to the notation on p.2.12 and 2.15 of
  Peakall and Smouse (2012).
  \item{phipt}{value of phiPT.}
  \item{vap}{variance among (between) populations (communities).}
  \item{n0}{standardisation factor N0, see  p.2.12 of Peakall and Smouse (2012).}
  \item{sst}{total distances sum of squares.}
  \item{ssg}{vector with two non-\code{NA} entriesm, within community
    sums of squares for communities \code{i} and \code{j}.}
  \item{msa}{mean square between communities.}
  \item{msw}{mean square within communities.}
}

    
\references{
  Bowcock, A. M., Ruiz-Linares, A., Tomfohrde, J., Minch, E., Kidd, J. R.,
  Cavalli-Sforza, L. L. (1994) High resolution of human evolutionary
  trees with polymorphic microsatellites. \emph{Nature} 368, 455-457.
  
  Cavalli-Sforza, L. L. and Edwards, A. W. F. (1967) Phylogenetic
  Analysis - Models and Estimation Procedures. \emph{The American Journal of
  Human Genetics} 19, 233-257.

  Meirmans, P. G. (2006) Using the AMOVA framework to estimate a
  standardized genetic differentiation measure. \emph{Evolution} 60, 2399-2402.
  
  Peakall, R. and Smouse P.E. (2012) GenAlEx Tutorial 2.
  \url{https://biology-assets.anu.edu.au/GenAlEx/Tutorials.html}
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{
  \code{\link{communitydist}}
}

\examples{
  options(digits=4)
  data(tetragonula)
  tnb <-
  coord2dist(coordmatrix=tetragonula.coord[83:120,],cut=50,file.format="decimal2",neighbors=TRUE)
  ta <- alleleconvert(strmatrix=tetragonula[83:120,])
  tai <- alleleinit(allelematrix=ta,neighborhood=tnb$nblist)
  tetracoms <-
  c(rep(1:3,each=3),4,5,rep(6:11,each=2),12,rep(13:19,each=2))
  phipt(tai,tetracoms,4,6)
  tdip <- diploidcomlist(tai,tetracoms,diploid=TRUE)
  cfchord(tdip[[4]],tdip[[6]])
  shared.problist(tdip[[4]],tdip[[6]])
  
}
\keyword{spatial}% __ONLY ONE__ keyword per line
\keyword{multivariate}% __ONLY ONE__ keyword per line
