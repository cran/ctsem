% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctACF.R
\name{ctResiduals}
\alias{ctResiduals}
\title{Extract Standardized Residuals from a ctsem Fit}
\usage{
ctResiduals(fit)
}
\arguments{
\item{fit}{A fitted model object generated by the ctsem package.}
}
\value{
A data table containing the standardized residuals for each subject and time point.
}
\description{
This function takes a fit object from the ctsem package and extracts the standardized residuals.
}
\details{
This function uses the \code{\link{ctStanKalman}} function to calculate the standardized residuals
and then extracts and formats them as a data table. The standardized residuals represent the differences
between the observed and predicted values, divided by the standard errors of the observations.
}
\examples{
data.table::setDTthreads(1) #ignore this line
# Example usage:
residuals <- ctResiduals(ctstantestfit)

}
\seealso{
\code{\link{ctStanKalman}}
}
