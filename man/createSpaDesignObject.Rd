% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createSpaDesignObject.R
\name{createSpaDesignObject}
\alias{createSpaDesignObject}
\title{Create spaDesign object}
\usage{
createSpaDesignObject(count_matrix, loc)
}
\arguments{
\item{count_matrix}{A gene expression count \code{matrix}}

\item{loc}{A \code{data.frame} containing the spatial coordinates and domain information for spots, corresponding to columns named as 'x', 'y', and 'domain'}
}
\value{
Returns a \code{spaDesign} object
}
\description{
Create spaDesign object
}
\examples{

## Create a spaDesign object
toyDATA  <- createSpaDesignObject(count_matrix = toyData$toyCount, loc = toyData$loc)

}
