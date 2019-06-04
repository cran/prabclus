\name{alleledist}
\alias{alleledist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Shared allele distance for diploid loci}
\description{
  Shared allele distance for codominant markers (Bowcock et al., 1994).
  One minus proportion of alleles shared by two individuals averaged
  over loci (loci with missing values for at least one individual are
  ignored).
}
\usage{
  alleledist(allelelist,ni,np,count=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{allelelist}{a list of lists. In the "outer" list, there are
    \code{np} lists, one for each locus. In the "inner" list, for every
    individual there is a vector of two codes (typically characters, see
    \code{\link{alleleinit}}) for the two alleles in that locus. Such a
    list can be constructed by \code{\link{unbuild.charmatrix}} out of
    the \code{charmatrix} component of an output object of
    \code{\link{alleleinit}}.}
  \item{ni}{integer. Number of individuals.}
  \item{np}{integer. Number of loci.}
  \item{count}{logical. If \code{TRUE}, the number of the individual to
    be processed is printed.}
}

\value{
  A symmetrical matrix of shared allele distances between individuals.
}

\references{
  Bowcock, A. M., Ruiz-Linares, A., Tomfohrde, J., Minch, E., Kidd, J. R.,
  Cavalli-Sforza, L. L. (1994) High resolution of human evolutionary
  trees with polymorphic microsatellites. \emph{Nature} 368, 455-457.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{
%-  \code{\link{chorddist}},\code{\link{neidist}},
  \code{\link{alleleinit}}, \code{\link{unbuild.charmatrix}}
}

\examples{
  data(tetragonula)
  tnb <-
  coord2dist(coordmatrix=tetragonula.coord[1:50,],cut=50,file.format="decimal2",neighbors=TRUE)
  ta <- alleleconvert(strmatrix=tetragonula[1:50,])
  tai <- alleleinit(allelematrix=ta,neighborhood=tnb$nblist,distance="none")
  str(alleledist((unbuild.charmatrix(tai$charmatrix,50,13)),50,13)) 
}
\keyword{cluster}% at least one, from doc/KEYWORDS
