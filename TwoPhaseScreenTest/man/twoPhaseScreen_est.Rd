
\name{twoPhaseScreen_est}
\alias{twoPhaseScreen_est}
\title{twoPhaseScreen_est}

\description{
To estimate the effect size of the screening test on the outcome.
}

\usage{
twoPhaseScreen_est(design.obj, data, interac.ind=c(0,0), brks=1.5) 
}
\arguments{
  \item{design.obj}{The object returned from the design function \code{\link{twoPhaseScreen_design}}.}
  \item{data}{A data frame containing two columns, where the first column are sample ids and 
              the second column binary (0-1) values for the screening test.}
  \item{interac.ind}{A binary vector of length two determining which models will have
                    the interaction term \code{x1*x2}. The four possibilities are: \cr
                    \code{c(1,1)} - added to both logistic and proportion hazards model \cr
                    \code{c(1,0)} - added to the logistic model only \cr
                    \code{c(0,1)} - added to the proportion hazards model only \cr
                    \code{c(0,0)} - not added to either model \cr
  }
  \item{brks}{Vector of cut-points for the piece-wise constant hazard function.}
}
\value{
  A list containing the vector of parameter estimates (\code{est}),
  the vector of standard errors (\code{se}), the covariance matrix (\code{cov}),
  \code{interac.ind}, \code{brks}, 
  and the results from \code{\link[stats]{optim}} (\code{optim.obj}).
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
  est.obj$est
  est.obj$se 
  }

}
