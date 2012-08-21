\name{prabclus-package}
\alias{prabclus-package}
%- Also NEED an `\alias' for EACH other topic documented here.
\docType{package}
\title{prabclus package overview}
\description{
  Here is a list of the main functions in package prabclus. Most other
  functions are auxiliary functions for these.
}
\section{Initialisation}{
  \describe{
    \item{prabinit}{Initialises presence/absence-, abundance- and
      multilocus data with dominant markers for use with most other key
      prabclus-functions.}
    \item{alleleinit}{Initialises multilocus data with codominant
      markers for use with key prabclus-functions.}
    \item{alleleconvert}{Generates the input format required by
      \code{\link{alleleinit}}.}
    }}
\section{Tests for clustering and nestedness}{
  \describe{
    \item{prabtest}{
      Computes the tests introduced in Hausdorf and Hennig (2003) and
      Hennig and Hausdorf (2004; these tests occur in some further 
      publications of ours but this one is the most detailed statistical
      reference) for presence/absence data. Allows use of the
      geco-dissimilarity (Hennig and Hausdorf, 2006).}
    \item{abundtest}{
      Computes the test introduced in Hausdorf and Hennig (2007) for
      abundance data.}
    \item{homogen.test}{A classical distance-based test for homogeneity
      going back to Erdos and Renyi (1960) and Ling (1973).}
    }}
\section{Clustering}{
  \describe{
    \item{prabclust}{Species clustering for biotic element analysis
      (Hausdorf and Hennig, 2007, Hennig and Hausdorf, 2004 and others),
      clustering of individuals for species delimitation (Hausdorf and
      Hennig, 2010) based on Gaussian mixture model clustering with
      noise as implemented in R-package \code{mclust}, Fraley and
      Raftery (1998),  on output of
      multidimensional scaling from distances as computed by
      \code{\link{prabinit}} or \code{\link{alleleinit}}. See also
      \code{\link{stressvals}} for help with choosing the number of
      MDS-dimensions.} 
    \item{hprabclust}{An unpublished alternative to
      \code{\link{prabclust}} using hierarchical clustering methods.}
    \item{lociplots}{Visualisation of clusters of genetic markers vs.
      clusters of species.}
    \item{NNclean}{Nearest neighbor based classification of observations
      as noise/outliers according to Byers and Raftery (1998).}
}}
\section{Dissimilarity matrices}{
  \describe{
    \item{alleledist}{Shared allele distance (see the corresponding help
      pages for references).}
    \item{dicedist}{Dice distance.}
    \item{geco}{geco coefficient, taking geographical distance into account.}
    \item{jaccard}{Jaccard distance.}
    \item{kulczynski}{Kulczynski dissimilarity.}
    \item{qkulczynski}{Quantitative Kulczynski dissimilarity for
      abundance data.}
}}
\section{Small conversion functions}{
  \describe{
    \item{coord2dist}{Computes geographical distances from geographical
      coordinates.}
    \item{geo2neighbor}{Computes a neighborhood list from geographical
      distances.}
    \item{alleleconvert}{A somewhat restricted function for conversion
      of different file formats used for genetic data with codominant markers.}
}}
\section{Data sets}{
      \code{\link{kykladspecreg}}, \code{\link{siskiyou}},
  \code{\link{veronica}}, \code{\link{tetragonula}}.
}

\references{
  Byers, S. and Raftery, A. E. (1998) Nearest-Neighbor Clutter
Removal for Estimating Features in Spatial Point Processes,
\emph{Journal of the American Statistical Association}, 93, 577-584.

  Erdos, P. and Renyi, A. (1960) On the evolution of random
  graphs. \emph{Publications of the Mathematical Institute of the
  Hungarian Academy of Sciences} 5, 17-61.

  Fraley, C.  and Raftery, A. E. (1998) How many clusters? Which
  clusterin method? - Answers via Model-Based Cluster Analysis.
  \emph{Computer Journal} 41, 578-588.

  Hausdorf, B. and Hennig, C. (2003) Nestedness of north-west European
land snail ranges as a consequence of differential immigration from
Pleistocene glacial refuges. \emph{Oecologia} 135, 102-109.

  Hausdorf, B. and Hennig, C. (2007) Null model tests of clustering of
  species, negative co-occurrence patterns and nestedness in meta-communities.
  \emph{Oikos} 116, 818-828. 

  Hausdorf, B. and Hennig, C. (2010) Species Delimitation Using Dominant
  and Codominant Multilocus Markers. \emph{Systematic Biology}, 59, 491-503.
  
  Hennig, C. and Hausdorf, B. (2004) Distance-based parametric bootstrap
tests for clustering of species ranges. \emph{Computational Statistics
  and
  Data Analysis} 45, 875-896.

  Hennig, C. and Hausdorf, B. (2006)  A robust distance coefficient
  between distribution areas incorporating geographic distances.
  \emph{Systematic Biology} 55, 170-175.

  Ling, R. F. (1973) A probability theory of cluster
  analysis. \emph{Journal of the American Statistical Association} 68,
  159-164.

}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche/}
}



