\name{alleleconvert}
\alias{alleleconvert}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Format conversion for codominant marker data}
\description{
  Codominant marker data (which here means: data with several diploid
  loci; two alleles per locus) can be represented in various ways. This
  function converts the formats \code{"genepop"} and \code{"structure"} into
  \code{"structurama"} and \code{"prabclus"}. \code{"genepop"}
  is a version of the format
  used by the package GENEPOP (Rousset, 2008), \code{"structure"} is a version
  of what is used by STRUCTURE (Pritchard et al., 2000), another one is
  \code{"structureb"}. \code{"structurama"}
  is a version of what is used by STRUCTURAMA (Huelsenbeck and
  Andolfatto, 2007) and \code{"prabclus"} is required by the function
  \code{\link{alleleinit}} in the present package.
}
\usage{
  alleleconvert(file=NULL,strmatrix=NULL, format.in="genepop",
                          format.out="prabclus",
                          alength=3,orig.nachar="000",new.nachar="-",
                          rows.are.individuals=TRUE, firstcolname=FALSE,
                          aletters=intToUtf8(c(65:90,97:122),multiple=TRUE),
                          outfile=NULL,skip=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{file}{string. Filename of input file, see details. One of
    \code{file} and \code{strmatrix} needs to be specified.}
  \item{strmatrix}{matrix or data frame of strings, see details. One of
    \code{file} and \code{strmatrix} needs to be specified.}
  \item{format.in}{string. One of \code{"genepop"},
    \code{"structure"}, or \code{"structureb"}, see details.}
  \item{format.out}{string. One of \code{"structurama"} or
    \code{"prabclus"}, see details.}
  \item{alength}{integer. If \code{format.in="genepop"}, length of code
    for a single allele.}
  \item{orig.nachar}{string. Code for missing values in input data.}
  \item{new.nachar}{string. Code for missing values in output data.}
  \item{rows.are.individuals}{logical. If \code{TRUE}, rows are
    interpreted as individuals and columns (variables if
    \code{strmatrix} is a data frame) as loci.}
  \item{firstcolname}{logical. If \code{TRUE}, it is assumed that the
    first column contains row names.}
  \item{aletters}{character vector. String of default characters for
    alleles if \code{format.out=="prabclus"} (the default is fine unless
    there is a locus that can have more than 62 different alleles in the
    dataset).}
  \item{outfile}{string. If specified, the output matrix (omitting
    quotes) is written to a file of this name (including row names if
    \code{fistcolname==TRUE}).}
  \item{skip}{number of rows to be skipped when reading data from a
      file (\code{skip}-argument of \code{\link{read.table}}).}
}
\details{
  The formats are as follows (described is the format within R, i.e.,
  for the input, the format of \code{strmatrix}; if \code{file} is
  specified, the file is read with
  \code{read.table(file,colClasses="character")} and should give the
  format explained below - note that \code{colClasses="character"}
  implies that quotes are not needed in the input file):
  \describe{
    \item{genepop}{Alleles are coded by strings of length \code{alength}
      and there is no space between the two alleles in a locus, so a
      value of \code{"258260"} means that in the corresponding locus the two
      alleles have codes 258 and 260.}
    \item{structure}{Alleles are coded by strings of arbitrary length.
      Two rows correspond to each inidividual, the first row containing
      the first alleles in all loci and the second row containing the
      second ones.}
    \item{structureb}{Alleles are coded by strings of arbitrary length.
      One row corresponds to each inidividual, containing first and
      second alleles in all loci (first and second allele of first locus, first
      and second allele of second locus etc.). This starts in the third
      row (first two have locus names and other information).}
    \item{structurama}{Alleles are coded by strings of arbitrary
      length. the two alleles in each locus are written with brackets
      around them and a comma in between, so \code{"258260"} in
      \code{"genepop"}
      corresponds to \code{"(258,260)"} in \code{"structurama"}.}
    \item{prabclus}{Alleles are coded by a single character and there is
      no space between the two alleles in a locus (e.g.,
      \code{"AC"}).}
  }
}

\value{
  A matrix of strings in the format specified as \code{format.out} with
  an attribute \code{"alevels"}, a vector of all used allele codes if
  \code{format.out=="prabclus"}, otherwise vector of allele codes of
  last locus.  
}
\references{
  Huelsenbeck, J. P., and P. Andolfatto (2007) Inference of population
  structure under a Dirichlet process model. \emph{Genetics} 175,
  1787-1802.
  
  Pritchard, J. K., M. Stephens, and P. Donnelly (2000) Inference of
  population structure using multi-locus genotype data. \emph{Genetics}
  155, 945-959.
  
  Rousset, F. (2008) genepop'007: a complete re-implementation of the
  genepop software for Windows and Linux. \emph{Molecular Ecology
  Resources} 8, 103-106.


}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{
  \code{\link{alleleinit}}
}

\examples{
  data(tetragonula)
# This uses example data file Heterotrigona_indoFO.dat
  str(alleleconvert(strmatrix=tetragonula))
  strucmatrix <-
    cbind(c("I1","I1","I2","I2","I3","I3"),
    c("122","144","122","122","144","144"),c("0","0","21","33","35","44"))
  alleleconvert(strmatrix=strucmatrix,format.in="structure",
    format.out="prabclus",orig.nachar="0",firstcolname=TRUE)
  alleleconvert(strmatrix=strucmatrix,format.in="structure",
    format.out="structurama",orig.nachar="0",new.nachar="-9",firstcolname=TRUE)
}
\keyword{manip}% __ONLY ONE__ keyword per line
