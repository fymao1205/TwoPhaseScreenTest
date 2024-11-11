
\name{twoPhaseScreen_design}
\alias{twoPhaseScreen_design}
\title{twoPhaseScreen_design}

\description{
To identify individuals to obtain screening test results from. This may be from testing stored biospecimens collected at a baseline visit.
}

\usage{
twoPhaseScreen_design(data, id, time.lower, time.upper, disease.at.baseline,
                      disease.after.baseline, covariate, phaseI.strata=NULL,
                      pcuts=c(1/3, 2/3), phaseII.design="RSD12", phaseII.nsamp=NULL, 
                      phaseII.weight=0.5) 
}
\arguments{
  \item{data}{Data frame containing all variables needed for analysis.}
  \item{id}{Column name in data containing unique patient ids.}
  \item{time.lower}{Column name in data referring to time vector defined as followed: For events, the time from baseline to the last disease-free ascertainment (with 0 used to denote prevalent events, no disease-free ascertainment, or the last disease-free ascertainment occurring at the baseline visit); for non-events, the time from baseline to the second-to-last disease-free ascertainment time (with 0 used to denote 0 or 1 disease-free ascertainments).}
  \item{time.upper}{Column name in data referring to another time vector defined as followed: For events, the time from baseline to the first time the event is detected (with 0 used to denote prevalent events); For non-events, this time from baseline to the last disease-free ascertainment (setting to 0 if there is no disease-free ascertainment or if the last disease-free ascertainment occurs at the baseline visit).}
  \item{disease.at.baseline}{Column name in data referring to disease status at baseline with coding as followed: 0 for disease-free, 1 for disease-detected, and 2 for unknown disease status.}
  \item{disease.after.baseline}{Column name in data referring to disease status after baseline with coding as followed: 0 for disease-free, 1 for disease-detected, and 2 for unknown disease status.}
  \item{covariate}{Column name in data referring to a binary (inexpensive) covariate. The referred column vector is coded as 0-1.}
  \item{phaseI.strata}{A character vector of variable names to define the strata. The
                      variable names must be a subset of \code{time.upper}, 
                      \code{disease.at.baseline}, \code{disease.after.baseline}, and 
                      \code{covariate}. The default is \code{disease.after.baseline}.}
  \item{pcuts}{An optional variable specifying a vector of cutpoints to discretize the continuous variable \code{time.upper}. 
              The default is \code{c(1/3, 2/3)}.}
  \item{phaseII.design}{Specified design with the following options: "RSD1" (univariate residual-dependent sampling design for \eqn{\beta_2}{beta_2}, the association parameter of the "expensive" screening test in the proportional hazards model), 
                                                                      "RSD2" (univariate residual-dependent sampling design for \eqn{\gamma_2}{gamma_2}, the association parameter of the "expensive" screening test in the logistic model), 
                                                                      "RSD12" (two-wave bivariate residual-dependent sampling design for \eqn{\beta_2}{beta_2} and \eqn{\gamma_2}{gamma_2}), 
                                                                      "srs" (simple random sampling design), 
                                                                      and "bal" (balanced-stratified sampling design).
	 
    If "RSD12" is used, we need to specify phaseII.weight to determine the proportion of the phase II sub-sample to be selected in wave 1 such that 1-phaseII.weight is selected in wave 2; The default is "RSD12".}
  \item{phaseII.nsamp}{Number of samples to select in phase II to measure
                       the "expensive" screening test.}
  \item{phaseII.weight}{Weight to be used in wave 1 selection with
                        \code{phaseII.design = "RSD12"}. The default is 0.5.}
}
\value{
  A list containing the objects \code{selected.ids}, \code{phaseI.strata.freq},
  \code{phaseII.strata.freq}, \code{data}, \code{options}, and \code{variable.map}.

  \code{selected.ids} is the vector of sample ids to collect the expensive screening test results from.

  \code{phaseI.strata.freq} is a frequency table reporting the size of each phaseI stratum.

  \code{phaseII.strata.freq} is a frequency table reporting the numbers selected in phase II from each phase I stratum.

  \code{data} is the data frame containing the phase I information including outcome-related (i.e. event status and associated time variables) variables and inexpensive covariates. 

  \code{options} is a list of options used in the analysis.

  \code{variable.map} is a named vector of the original and transformed variable names. 

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

ret <- twoPhaseScreen_design(data_interval, id, time.lower, time.upper, disease.at.baseline,
                         disease.after.baseline, covariate, phaseI.strata=c("d1", "d2", "x1"),
                         pcuts=c(1/3, 2/3), phaseII.design="RSD12", phaseII.nsamp=2000, 
                         phaseII.weight=0.5)
ret$phaseI.strata.freq
ret$phaseII.strata.freq
length(ret$selected.ids)

}
