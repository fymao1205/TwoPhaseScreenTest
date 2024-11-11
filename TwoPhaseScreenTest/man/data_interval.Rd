
\name{data_interval}
\alias{data_interval}
\alias{data_interval}

\title{
  Data for examples
}
\description{
 A simulated data example for left- and interval-censored data.
}
\details{
  A data frame containing a simulated complete data with interval-censored outcomes. It contains \cr
                    \code{id} - patient id \cr
                    \code{z} - true status of prevalent event \cr
                    \code{t} - true event time (positive values if incident events; 0 if prevalent event) \cr
                    \code{r.z} - indicator of observing z at baseline \cr
                    \code{z.o} - observed values of \eqn{z}{z}; z.o=2 if missing \cr
                    \code{a1} - For events, the time from baseline to the last disease-free ascertainment (with 0 used to denote prevalent events, no disease-free ascertainment, or the last disease-free ascertainment occurring at the baseline visit); for non-events, the time from baseline to the second-to-last disease-free ascertainment time (with 0 used to denote 0 or 1 disease-free ascertainments) \cr
                    \code{a2} - For events, the time from baseline to the first time the event is detected (with 0 used to denote prevalent events); For non-events, the time from baseline to the last disease-free ascertainment (setting to 0 if there is no disease-free ascertainment, or the last disease-free ascertainment occurring at the baseline visit) \cr
                    \code{u} - lower bound of the observed time interval that brackets the 
                             latent event time \eqn{t}{t} \cr
                    \code{v} - upper bound of the observed time interval that brackets the 
                      latent event time \eqn{t}{t} \cr
                    \code{d1} - observed disease status at baseline; 0 for disease-free, 1 for disease-detected, and 2 for unknown disease status \cr
                    \code{d2} - observed disease status after baseline; 0 for disease-free, 1 for disease-detected, and 2 for unknown disease status \cr
                    \code{x1} - binary "inexpensive" covariate \cr
                    \code{x2} - binary "expensive" covariate \cr
}
\examples{

 data(data_interval, package="TwoPhaseScreenTest")

 # Display a few rows of the data
 data_interval[1:5, ]

}
\keyword{data}
%\seealso{ \code{\link{twoPhaseScreenTest}} }
