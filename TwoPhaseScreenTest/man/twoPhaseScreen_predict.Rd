
\name{twoPhaseScreen_predict}
\alias{twoPhaseScreen_predict}
\title{twoPhaseScreen_predict}

\description{
Compute the positive predicted value (PPV(t)) and the complement of the negative
predicted value (cNPV(t)) of a screening test for cumulative outcomes as of a specified time point.
}

\usage{
twoPhaseScreen_predict(est.obj, t=0, grad=1e-06) 
}
\arguments{
  \item{est.obj}{The object returned from the estimation function \code{\link{twoPhaseScreen_est}}.}
  \item{t}{Time since screening test to calculate PPV(t) and cNPV(t). The default is 0.}
  \item{grad}{A small positive value for computing numerical derivatives. The default is 1e-6.}
}
\value{
  A list containing the estimated cumulative risks and corresponding standard errors for screening test
  values of 0 (\code{est.cNPV, se.cNPV}) and 1 (\code{est.PPV, se.PPV}).
}

\details{
  See the manuscript for complete details of the methods.
}

\examples{
  data(data_interval, package="TwoPhaseScreenTest")

  id                     <- "id"
  time.lower             <- "a1"
  time.upper             <- "a2"
  disease.at.baseline    <- "d1"
  disease.after.baseline <- "d2"
  covariate              <- "x1"

  design.obj <- twoPhaseScreen_design(data_interval, id, time.lower, time.upper, disease.at.baseline,
                         disease.after.baseline, covariate, phaseI.strata=c("d1", "d2", "x1"),
                         pcuts=c(1/3, 2/3), phaseII.design="RSD12", phaseII.nsamp=2000, 
                         phaseII.weight=0.5)

  expensive_cov_data <- data_interval[, c("id", "x2")]
  \donttest{
  est.obj <- twoPhaseScreen_est(design.obj, expensive_cov_data, interac.ind=c(0,0), brks=1.5) 
  twoPhaseScreen_predict(est.obj, t=0)
  twoPhaseScreen_predict(est.obj, t=3)
  }
}
