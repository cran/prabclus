\name{tetragonula}
\alias{tetragonula}
\alias{tetragonula.coord}
\docType{data}
% \non_function{}
\title{Microsatellite genetic data of Tetragonula bees}
\description{
  Genetic data for 236 Tetragonula (Apidae) bees from Australia
  and Southeast Asia, see  Franck et al. (2004). The data give pairs of
  alleles (codominant markers) for 13 microsatellite loci.
}
\usage{data(tetragonula)}
\format{
  Two objects are generated:
  \describe{
    \item{tetragonula}{A data frame with 236 observations and 13 string
 variables. Strings consist of six digits each. The
 format is derived from the data format used by the software GENEPOP
 (Rousset 2010). Alleles have a three digit code, so a value of
 \code{"258260"} on variable V10 means that on locus 10 the two alleles have
 codes 258 and 260. \code{"000"} refers to missing values.}
\item{tetragonula.coord}{a 236*2 matrix. Coordinates of locations of
    individuals in decimal format, i.e. the first number is latitude
    (negative values are South), with minutes and seconds
    converted to fractions. The second number is longitude (negative
    values are West).}
  }
}
\source{
  Franck, P., E. Cameron, G. Good, J.-Y. Rasplus, and B. P. Oldroyd
  (2004) Nest architecture and genetic differentiation in a species
  complex of Australian stingless bees. \emph{Mol. Ecol.} 13, 2317-2331.

  Rousset, F. (2010) Genepop 4.0 for Windows and Linux.
  \url{http://kimura.univ-montp2.fr/~rousset/Genepop.pdf}
}
\details{
   Reads from example data file \code{Heterotrigona_indoFO.dat}. 
   }

\examples{
data(tetragonula)
}
\keyword{datasets}
