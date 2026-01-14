\name{veronica}
\alias{veronica}
\alias{veronica.coord}
\docType{data}
% \non_function{}
\title{Genetic AFLP data of Veronica plants}
\description{
  0-1 data indicating whether dominant markers are present for 583
  different AFLP bands ranging from 61 to 454 bp of 207 plant
  individuals of Veronica (Pentasepalae) from the Iberian Peninsula and
  Morocco (Martinez-Ortega et al., 2004). 
}
\usage{data(veronica)}
\format{
  Two objects are generated:
  \describe{
  \item{veronica}{0-1 matrix with 207 individuals (rows) and 583 AFLP
    bands (columns).}
  \item{veronica.coord}{a 207*2 matrix. Coordinates of locations of
    individuals in decimal format, i.e. the first number is latitude
    (negative values are South), with minutes and seconds
    converted to fractions. The second number is longitude (negative
    values are West).}
  }
}
\source{
   Martinez-Ortega, M. M., L. Delgado, D. C. Albach,
   J. A. Elena-Rossello, and E. Rico (2004). Species boundaries
   and phylogeographic patterns in cryptic taxa inferred from AFLP
   markers: Veronica subgen. Pentasepalae (Scrophulariaceae) in the
   Western Mediterranean.\emph{Syst. Bot.} 29, 965-986.
}
\details{
  Reads from example data files \code{MartinezOrtega04AFLP.dat,
    MartinezKoord.dat}. 
  }

\examples{
data(veronica)
}
\keyword{datasets}
