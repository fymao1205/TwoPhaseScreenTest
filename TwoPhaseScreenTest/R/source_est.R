
# using design...
# collect x2 for the phase 2 sub-sample ...
# now for data analysis ...

twoPhaseScreen_est <- function(design.obj, data, interac.ind=c(0,0), brks=1.5) {

  # Check arguments
  check_design.obj(design.obj)
  check_data_est(data)
  check_interac.ind(interac.ind) 
  check_num(brks, "brks", len=0, pos=1)  

  objlist <- list(interac.ind=interac.ind, brks=brks)

  ret <- tps_est_main(design.obj, data, objlist)  

  ret
}

tps_est_main <- function(design.obj, data, objlist) {

  dt_ph2 <- design.obj$data
  idvar  <- "id" # Renamed to id in design function, id gives row numbers
  if (!(idvar %in% colnames(dt_ph2))) {
    msg <- paste0("ERROR: ", idvar, " not found in design.obj$data")
    stop(msg)
  }

  # Initialize expensive covariate
  dt_ph2[, "x2"] <- NA

  # Merge in expensive covariate
  rnms <- rownames(dt_ph2)
  rows <- match(trimws(rnms), trimws(data[, 1, drop=TRUE]))
  tmp  <- !is.na(rows)
  rows <- rows[tmp]
  if (!length(rows)) stop("ERROR: no sample ids match in data")
  dt_ph2[tmp, "x2"] <- as.numeric(data[rows, 2, drop=TRUE])

  tmp               <- !(dt_ph2[, idvar, drop=TRUE] %in% design.obj$selected.ids)
  dt_ph2[tmp, "x2"] <- NA
  dt_ph2[, "r"]     <- 1
  dt_ph2[tmp, "r"]  <- 0

  # Check for error, expensive cov must be binary for selected subjects
  tmp2 <- !(dt_ph2[!tmp, "x2", drop=TRUE] %in% 0:1)
  if (any(tmp2)) stop("ERROR: selected.ids must contain binary data for the expensive covariate")
  rm(tmp, tmp2, rows); gc()

  # if x1*x2 interaction term was added to the two models
  # (1,1): added to both
  # (1,0): added to logistic but PH
  # (0,1): added to the PH but logistic
  # (0,0): none
  interac.ind <- objlist$interac.ind
 
  # cut-points for the piece-wise constant hazard function
  brks <- objlist$brks

  # estimation based on maximum likelihood
  est <- est.obs_pwc.f(dt_ph2, interac.ind, brks)

  # inference: variance matrix 
  ase <- ase_pwc.f(grad=1e-06, est, dt_ph2, interac.ind, brks)

  # Get names for se, cov
  nms   <- tps_getCovNames(interac.ind, brks)
  parms <- est$par
  se    <- ase$ase
  cov   <- ase$avar
  if (length(nms) == length(se)) {
    names(parms)  <- nms
    names(se)     <- nms
    rownames(cov) <- nms
    colnames(cov) <- nms
  }

  # Include interac.ind and brks in return object
  list(est=parms, se=se, cov=cov, interac.ind=interac.ind, brks=brks, optim.obj=est)
}

tps_getCovNames <- function(interac.ind, brks) {

  if(any(is.finite(brks))){
    len = length(brks)+1
  }else{
    len = 1
  }
  
  if(all(interac.ind==c(0,0))){
    len.beta=2
    len.gam=3
  }else if(all(interac.ind==c(1,0))){
    len.beta=3
    len.gam=3
  }else if(all(interac.ind==c(0,1))){
    len.beta=2
    len.gam=4
  }else{
    len.beta=3
    len.gam=4
  }
  
  #eta = para[1:2]
  #p1 = expit.f(para[3])
  #alp = exp(para[3+1:len])
  #beta = c(para[3+len+1:len.beta], 0)
  #gam = c(para[3+len+len.beta+1:len.gam],0)
  
  ret <- c("eta1", "eta2", "p1", paste0("alpha", 1:len), 
           paste0("beta", 1:len.beta), paste0("gamma", 1:len.gam))
  ret
}

