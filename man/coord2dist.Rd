\name{coord2dist}
\alias{coord2dist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Geographical coordinates to distances}
\description{
  Computes geographical distances from geographical coordinates
}
\usage{
  coord2dist(file=NULL, coordmatrix=NULL, cut=NULL,
                       file.format="degminsec",
                       output.dist=FALSE, radius=6378.137,
              fp=1/298.257223563, neighbors=FALSE)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{file}{string. A filename for the coordinate file. The file
    should have 2, 4 or 6 numeric columns and one row for each location.
    See \code{file.format}. One of \code{file} and \code{coordmatrix}
    needs to be specified (if \code{coordmatrix} is not specified,
    coordinates are read from \code{file}).}
  \item{coordmatrix}{something that can be coerced into a matrix with
    2, 4 or 6 columns. Matrix of coordinates, one row for each
    location. See \code{file.format}. One of \code{file} and
    \code{coordmatrix} needs to be specified.}
  \item{cut}{numeric. Only active if \code{neighbors==TRUE}; see
    \code{neighbors}.}
  \item{file.format}{one of \code{"degminsec"}, \code{"decimal2"} or
    \code{"decimal4"}. The format of the required file or
    \code{coordmatrix} consists of the following columns:
    \describe{
      \item{"degminsec"}{6 columns; the first three give degrees,
	minutes and seconds for latitude, columns 4-6 the same for
	longitude.  Values in column 1 and 4 can be positive or negative
	(negative means "South", "West", respectively). Values in the
	other columns should be non-negative.}
      \item{"decimal2"}{2 columns; the first one gives latitude, the
	second one longitude in proper decimal notation. Values can be
	positive or negative (negative means "South", "West",
	respectively).}
      \item{"decimal4"}{4 columns; the first two give latitude, no. 3
	and 4 give longitude. Values in column 1 and 3 can be positive
	or negative (negative means "South", "West",
	respectively). The give integer degrees. Values in the other
	columns should be non-negative. They give percentages
	(\code{<=100}).}
  }}
  \item{output.dist}{logical. If \code{TRUE}, the resulting distance
    matrix is given out as a \code{\link{dist}} object.}
  \item{radius}{numeric. Radius of the earth in km used in computation
    (the default is the equatorial radius but this is not the uniquely
    possible choice).}
  \item{fp}{flattening of the earth; the default is from WGS-84.}
  \item{neighbors}{logical. If \code{TRUE}, a neighborhood list is also
    computed, listing for every location all locations with distance
    \code{<=cut} as neighbors.}      
}
\value{
  If \code{neighbors==TRUE}, a
  list with components
  \item{distmatrix}{distance matrix between locations. See
    \code{output.dist} above. This is in km by default; the measurement
    unit is determined by the value used for \code{radius}.}
  \item{nblist}{list with a vector for every location containing the
    numbers of its neighbors, see \code{neighbors}.}
  If \code{neighbors==FALSE}, only the distance matrix.
}
\references{
  German Wikipedia from 29 August 2010:
  \url{http://de.wikipedia.org/wiki/Orthodrome}
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{
  \code{\link{geo2neighbor}}
}

\examples{
  options(digits=4)
  data(veronica)
  coord2dist(coordmatrix=veronica.coord[1:20,], cut=20, file.format="decimal2",neighbors=TRUE)
}
\keyword{math}% __ONLY ONE__ keyword per line
