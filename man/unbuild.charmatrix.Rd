\name{unbuild.charmatrix}
\alias{unbuild.charmatrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Internal: create allele list out of character matrix}
\description{
  Creates a list of lists, such as required by \code{\link{alleledist}},
  from the \code{charmatrix} component of an
  \code{\link{alleleobject}}. 
}

\usage{
   unbuild.charmatrix(charmatrix,n.individuals,n.variables)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{charmatrix}{matrix of characters in which there are two rows for
    every individual corresponding to the two alleles in every locus
    (column). Entries are allele codes but missing values are coded as
    \code{NA}.}  
  \item{n.individuals}{integer. Number of individuals.}
  \item{n.variables}{integer. Number of loci.}
}
  
\value{
  A list of lists. In the "outer" list, there are
  \code{n.variables} lists, one for each locus. In the "inner" list, for every
    individual there is a vector of two codes (typically characters, see
    \code{\link{alleleinit}}) for the two alleles in that locus.
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{alleleinit}}, \code{\link{build.charmatrix}}}

\examples{
  data(tetragonula)
  tnb <-
  coord2dist(coordmatrix=tetragonula.coord[1:50,],cut=50,file.format="decimal2",neighbors=TRUE)
  ta <- alleleconvert(strmatrix=tetragonula[1:50,])
  tai <- alleleinit(allelematrix=ta,neighborhood=tnb$nblist,distance="none")
  str(unbuild.charmatrix(tai$charmatrix,50,13)) 
}
\keyword{manip}% at least one, from doc/KEYWORDS




