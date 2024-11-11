RtoC_logLn_cal1_wei <- function(d1, d2, a1, a2, x1, x2, gam, alp, beta) {

  d1_size   <- length(d1)
  gam_size  <- length(gam)
  alp_size  <- length(alp)
  beta_size <- length(beta)
  ret       <- rep(-9999, d1_size)

  # Check x2 and set it to a vector of ones if needed
  if (length(x2) != d1_size) x2 <- rep(1, d1_size)

  tmp <- .C("C_logLn_cal1_wei", as.integer(d1), as.integer(d2), 
            as.numeric(a1), as.numeric(a2), as.numeric(x1), as.numeric(x2), 
            as.numeric(gam), as.numeric(alp), as.numeric(beta),
            as.integer(d1_size), as.integer(gam_size), as.integer(alp_size), 
            as.integer(beta_size), ret=as.numeric(ret),
            NAOK=TRUE, PACKAGE="TwoPhaseScreenTest")
  ret <- tmp$ret

  ret
}

RtoC_logLn_obs_pwc <- function(r, d1, d2, a1, a2, x1, x2, gam, alp, beta, eta, p1, brks) {

  d1_size   <- length(d1)
  gam_size  <- length(gam)
  alp_size  <- length(alp)
  beta_size <- length(beta)
  eta_size  <- length(eta)
  brks_size <- length(brks)
  ret       <- rep(-9999, d1_size)

  tmp <- .C("C_logLn_obs_pwc", as.integer(r), as.integer(d1), as.integer(d2), 
            as.numeric(a1), as.numeric(a2), as.numeric(x1), as.numeric(x2), 
            as.numeric(gam), as.numeric(alp), as.numeric(beta), as.numeric(eta), 
            as.numeric(p1), as.numeric(brks),
            as.integer(d1_size), as.integer(gam_size), as.integer(alp_size), 
            as.integer(beta_size), as.integer(eta_size), as.integer(brks_size), 
            ret=as.numeric(ret), NAOK=TRUE, PACKAGE="TwoPhaseScreenTest")
  ret <- tmp$ret

  ret
}

RtoC_logL_cal_pwc <- function(d1, d2, t1, t2, x1, x2, gam, alp, beta, eta, p1, brks) {

  ngam  <- length(gam)
  nalp  <- length(alp)
  nbeta <- length(beta)
  neta  <- length(eta)
  nbrks <- length(brks)
  ret   <- -99999

  tmp <- .C("C_logL_cal_pwc", as.numeric(d1), as.numeric(d2), as.numeric(t1), as.numeric(t2), 
            as.numeric(x1), as.numeric(x2), as.numeric(gam), as.numeric(alp), as.numeric(beta),
            as.numeric(eta), as.numeric(p1),  as.numeric(brks), as.integer(ngam),
            as.integer(nalp), as.integer(nbeta), as.integer(neta), as.integer(nbrks), 
            ret=as.numeric(ret), PACKAGE="TwoPhaseScreenTest")
  ret <- tmp$ret

  ret
}

RtoC_ppwc <- function(q, cuts, levels, lower, logInd) {

  ncuts   <- length(cuts)
  nlevels <- length(levels)
  ret     <- -9999

  tmp <- .C("C_ppwc", as.numeric(q), as.numeric(cuts), as.numeric(levels), as.integer(lower),
            as.integer(logInd), as.integer(ncuts), as.integer(nlevels), 
            ret=as.numeric(ret), PACKAGE="TwoPhaseScreenTest")
  ret <- tmp$ret
  ret
}
