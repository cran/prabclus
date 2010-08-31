\name{lociplots}
\alias{lociplots}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Visualises clusters of markers vs. species}
\description{
  Given a clustering of individuals from \code{\link{prabclust}} (as
  generated in species delimitation) and a clustering of markers (for
  example dominant markers of genetic loci), \code{lociplots} visualises 
  the presence of markers against the clustering of individuals and
  computes some statistics.   
}
\usage{
  lociplots(indclust,locclust,locprab,lcluster,
                      symbols=NULL,brightest.grey=0.8,darkest.grey=0,
                      mdsdim=1:2)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{indclust}{\code{\link{prabclust}}-object. Clustering of individuals.}
  \item{locclust}{vector of integers. Clustering of markers/loci.}
  \item{locprab}{\code{\link{prab}}-object in which the markers are what
    the help page of \code{\link{prabinit}} refers to as "species"
    (i.e., reverse of what is used for species delimitation clustering;
    for data sets with codominant markers, such an object can be
    constructed by use of \code{\link{allele2zeroone}} before
    \code{\link{prabinit}}.)}
  \item{lcluster}{integer. Number of cluster in \code{locclust} for
    which plot and statistics are produced.}
  \item{symbols}{vector of plot symbols. If \code{NULL},
    \code{indclust$symbols} is used.}
  \item{brightest.grey}{numeric between 0 and 1. Brightest grey value
    used in plot for individuals with smallest marker percentage, see
    details.}
  \item{darkest.grey}{numeric between 0 and 1. Darkest grey value
    used in plot for individuals with highest marker percentage, see
    details.}
  \item{mdsdim}{vector of two integers. The two MDS variables taken from
    \code{indclust} used for visualisation.}  
}

\details{
  Plot and statistics are based on the individual marker percentage,
  which is the percentage of markers present in an individual of the
  markers belonging to cluster no. \code{lcluster}. In the plot, the
  grey value visualises the marker percentage.
}
  
\value{
  list with components
  \item{locfreq}{vector of individual marker percentages.}
  \item{locfreqmin}{vector of minimum individual marker precentages for
    each cluster in \code{indclust}-clustering (the first value refers
    to the "noise component", if present).}
  \item{locfreqmax}{vector of maximum individual marker precentages for
    each cluster in \code{indclust}-clustering (the first value refers
    to the "noise component", if present).}
  \item{locfreqmean}{vector of average individual marker precentages for
    each cluster in \code{indclust}-clustering (the first value refers
    to the "noise component", if present).}  
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{prabclust}}}

\examples{
  data(veronica)
  vei <- prabinit(prabmatrix=veronica[1:50,],distance="jaccard")
  ppv <- prabclust(vei)
  veloci <- prabinit(prabmatrix=veronica[1:50,],rows.are.species=FALSE)
  velociclust <- prabclust(veloci,nnk=0)
  lociplots(ppv,velociclust$clustering,veloci,lcluster=3)
}
\keyword{cluster}% at least one, from doc/KEYWORDS




