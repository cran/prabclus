\name{prabinit}
\alias{prabinit}
\alias{print.prab}
\alias{prab}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Presence-absence/abundance matrix initialization}
\description{
  \code{prabinit} converts a matrix into an object
  of class \code{prab} (presence-absence). The matrix may be read from a
  file or an R-object. It may be a 0-1 matrix or a matrix with
  non-negative entries (usually abundances).
  \code{print.prab} is a print method for such
  objects.

  Documentation here is in terms of biotic elements analysis (species
  are to be clustered). For species delimitation with dominant markers,
  see Hausdorf and Hennig (2010), individuals take the role of species
  and loci take the role of regions.
}
\usage{
prabinit(file = NULL, prabmatrix = NULL, rows.are.species = TRUE,
neighborhood = "none", nbbetweenregions=TRUE, geodist=NULL, gtf=0.1,
distance = "kulczynski", toprab = FALSE, toprabp
= 0.05, outc = 5.2)

\method{print}{prab}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{file}{string. non-negative matrix ASCII file (such as example dataset
    \code{kykladspecreg.dat}) from which the
    matrix is read by \code{read.table}. The usual interpretation is
    that it is a species-by-regions matrix of species presences/absences
    (0-1 matrix) or abundances.}
  \item{prabmatrix}{matrix with non-negative entries. Either \code{file}
    or \code{prabmatrix} should be \code{NA}.}
  \item{rows.are.species}{logical. If \code{TRUE}, rows are interpreted
    as species and columns are interpreted as regions. In this case,
    rows and columns are interchanged by \code{prabinit}.}
  \item{neighborhood}{A string or a list with a component for
    every region. The
    components are vectors of integers indicating
    neighboring regions. A region without neighbors (e.g., an island)
    should be assigned a vector \code{numeric(0)}. If \code{neighborhood}
    is a filename, it is attempted to read such a list from a file,
    where every row should correspond to one region (such as example
    dataset \code{nb.dat}). If
    \code{neighborhood="none"}, all neighborhoods are set to
    \code{numeric(0)}. The neighborhood can be tested by
    \code{\link{nbtest}} for consistency.}
  \item{nbbetweenregions}{logical. If \code{TRUE}, the neighborhood is
    defined between regions as explained above. Otherwise it is defined
    between species (or individuals, if this is used for species
    delimitation).}
  \item{geodist}{matrix of non-negative reals. Geographical distances
    between regions. Only used if \code{distance="geco"}.}
  \item{gtf}{tuning constant for geco-distance if \code{distance="geco"},
    see \code{geco}.}
  \item{distance}{\code{"kulczynski"}, \code{"jaccard"}, \code{"geco"},
    \code{"qkulczynski"}, \code{"logkulczynski"} (this calls function
    \code{qkulczynski} with \code{log.distance=TRUE}), \code{"dice"}, or
    \code{"none"}. The distance measure
    between species to compute by \code{prabinit}.}
  \item{toprab}{logical. If \code{TRUE}, a presence-absence matrix is
    computed from the non-negative input matrix. "Absence", i.e., the
    entry 0, is chosen if the original entry is 0, or the original entry
    is smaller than or equal to \code{toprabp} times the
    sum of entries in the
    corresponding region, and log(original entry) is considered to be a
    lower outlier compared with the other entries of the corresponding
    species (see \code{outc}). "Presence", i.e., the entry 1, thus
    means that the original
    entry is non-negligible w.r.t. the species or w.r.t. the region.}
  \item{toprabp}{numerical between 0 and 1, see \code{toprab}.}
  \item{outc}{numerical. Tuning constant for the outlier identification
    associated with \code{toprab=TRUE}. An entry smaller than or equal to
    \code{outc*mad} times the median is considered as a lower outlier.}
  \item{x}{object of class \code{prab}.}
  \item{...}{necessary for print method.}
}

\details{
  Species that are absent in all regions are omitted.}
 
\value{
  \code{prabinit} produces
  an object of class \code{prab}, which is a list with components
  \item{distmat}{distance matrix between species.}
  \item{prab}{abundance or presence/absence matrix (if presence/absence,
    the entries are logical). Rows are regions, columns are
    species.}
  \item{nb}{neighborhood list, see above.}
  \item{regperspec}{vector of the number of regions occupied by a
    species.}
  \item{specperreg}{vector of the number of species present in a region.}
  \item{n.species}{number of species (in the \code{prab}-object, see 
    \code{nonzero}).}
  \item{n.regions}{number of regions.}
  \item{distance}{string denoting the chosen distance measure.}
  \item{geodist}{non-negative matrix. see above.}
  \item{gtf}{numeric. see above.}
  \item{spatial}{\code{TRUE}, if there is a specified neighborhood
    structure.} 
  \item{nonempty.species}{logical vector. The length is the number of species
    in the original file/matrix. If \code{FALSE}, the corresponding species 
    had only zero entries and was therefore absent. Note that these species
    are not included in any other component of a \code{prab} object, i.e.,
    \code{n.species} is the number of \code{TRUE}-entries in
    \code{nonzero}.}
  \item{nbbetweenregions}{see above.}
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\references{
  Hausdorf, B. and Hennig, C. (2010) Species Delimitation Using Dominant
  and Codominant Multilocus Markers. \emph{Systematic Biology},
  published online 
  \url{http://sysbio.oxfordjournals.org/cgi/content/full/syq039?ijkey=LTDbmfF3LhbOOQY&keytype=ref}
}

\seealso{\code{\link{read.table}}, \code{\link{jaccard}},
  \code{\link{kulczynski}}, \code{\link{geco}},
  \code{\link{qkulczynski}}, \code{\link{nbtest}},
  \code{\link{alleleinit}}
  }

\examples{
# If you want to use your own ASCII data files, use
# x <- prabinit(file="path/prabmatrixfile",
# neighborhood="path/neighborhoodfile")
data(kykladspecreg)
data(nb)
prabinit(prabmatrix=kykladspecreg, neighborhood=nb)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
