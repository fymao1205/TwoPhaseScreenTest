
\name{data_csd}
\alias{data_csd}
\alias{data_csd}

\title{
  Data for examples
}
\description{
 A simulated data example for current status data.
}
\details{
  A data frame containing a simulated complete data with current-status observational schemes. It contains \cr
                    \code{id} - patient id \cr
                    \code{z} - true status of prevalent event \cr
                    \code{t} - true event time (positive values if incident events; 0 if prevalent event) \cr
                    \code{r.z} - indicator of observing z at baseline \cr
                    \code{a1} - a1=0 for all \cr
                    \code{a2} - For events, the time from baseline to the first time the event is detected (with 0 used to denote prevalent events); For non-events, the time from baseline to the last disease-free ascertainment (setting to 0 if there is no disease-free ascertainment, or the last disease-free ascertainment occurring at the baseline visit) \cr
                    \code{u} - lower bound of the observed time interval that brackets the 
                              latent event time \eqn{t}{t} \cr
                    \code{v} - upper bound of the observed time interval that brackets the
                              latent event time \eqn{t}{t} \cr
                    \code{d1} - observed disease status at baseline; 0 for disease-free, 1 for disease-detected, and 2 for unknown status \cr
                    \code{d2} - observed disease status after baseline; 0 for disease-free, 1 for disease-detected, and 2 for unknown status \cr
                    \code{x1} - binary "inexpensive" covariate \cr
                    \code{x2} - binary "expensive" covariate \cr
}
\examples{

 data(data_csd, package="TwoPhaseScreenTest")

 # Display a few rows of the data
 data_csd[1:5, ]

}
\keyword{data}
%\seealso{ \code{\link{twoPhaseScreenTest}} }
