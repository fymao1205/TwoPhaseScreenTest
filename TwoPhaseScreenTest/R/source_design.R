twoPhaseScreen_design <- function(data, id, time.lower, time.upper, disease.at.baseline,
                         disease.after.baseline, covariate, phaseI.strata=NULL,
                         pcuts=c(1/3, 2/3), phaseII.design="RSD12", phaseII.nsamp=NULL, 
                         phaseII.weight=0.5) {

  check_data(data)
  check_col(id, "id", data, num=0)
  check_col(time.lower, "time.lower", data, nonNeg=1)
  check_col(time.upper, "time.upper", data, nonNeg=1)
  check_col(disease.at.baseline, "disease.at.baseline", data, vals=0:2)
  check_col(disease.after.baseline, "disease.after.baseline", data, vals=0:2)
  check_col(covariate, "covariate", data, num=1, bin=1)
  pcuts <- check_cuts(pcuts)
  obj   <- c(disease.after.baseline, disease.at.baseline, time.upper, covariate)
  phaseI.strata <- check_strata.cols(phaseI.strata, data, obj)
  check_phaseII.design(phaseII.design)
  check_n2.samp(phaseII.nsamp, data)    
  check_num(phaseII.weight, "phaseII.weight", len=1, pos=1, min=0, max=1)  

  objlist <- list(id=id, time.lower=time.lower, time.upper=time.upper,  
              disease.at.baseline=disease.at.baseline, 
              disease.after.baseline=disease.after.baseline, covariate=covariate,
              pcuts=pcuts, phaseI.strata=phaseI.strata,
              phaseII.design=phaseII.design, phaseII.nsamp=phaseII.nsamp, 
              phaseII.weight=phaseII.weight)

  # Call main function
  ret <- tps_design_main(data, objlist)
  
  # Set return list
  ret <- tps_design_setRet(ret)

  ret
}

tps_design_setRet <- function(retlist) {

  # Return data, phaseI counts, phaseII counts, phaseII ids,
  #   variable map

  objlist <- retlist$objlist

  # PhaseI data
  data <- retlist$phaseI$dt_ext

  # Get phaseI freq counts
  svar  <- "group"
  freq1 <- table(data[, svar, drop=TRUE])

  # Rename group variable
  vars <- objlist$phaseI.strata
  new  <- paste0(vars, collapse="_")
  new  <- paste0("strata.", new)
  cx   <- colnames(data)
  tmp  <- cx %in% svar
  if (any(tmp)) {
    svar    <- new
    cx[tmp] <- new
  }
  colnames(data) <- cx

  # Get selected ids, the row numbers are returned in phaseII object
  ids2   <- retlist$phaseII$s_id
  rnms   <- rownames(data)
  selids <- rnms[ids2]

  # Get phaseII freq counts
  tmp   <- rnms %in% selids
  freq2 <- table(data[tmp, svar, drop=TRUE])

  # Variable map
  vmap <- objlist$map.origToNew

  # options
  nms <- c("pcuts", "phaseI.strata", "phaseII.design", "phaseII.nsamp",
           "phaseII.weight")
  op  <- objlist[nms]

  list(selected.ids=selids, phaseI.strata.freq=freq1, phaseII.strata.freq=freq2,
       data=data, options=op, variable.map=vmap)
}

tps_phaseII.design <- function(x) {

  ret <- x
  if (x == "RSD12") {
    ret <- "BRSD-eff-seq"
  } else if (x == "RSD1") {
    ret <- "RSD-Smu1"
  } else if (x == "RSD2") {
    ret <- "RSD-Smu2"
  }
  ret
}

tps_design_main <- function(data, objlist) {

  # Set up data
  tmp     <- tps_design_setData(data, objlist) 
  data    <- tmp$data
  objlist <- tmp$objlist
  rm(tmp); gc()

  # Get correct phaseII.design given the name changes
  objlist$phaseII.design.orig <- tps_phaseII.design(objlist$phaseII.design)

  # Set or check n2.samp since rows of data might have been removed
  tmp <- objlist[["phaseII.nsamp", exact=TRUE]]
  if (!length(tmp)) {
    objlist$phaseII.nsamp <- floor(nrow(data)/2)
  } else {
    check_n2.samp(tmp, data)
  }

  # Set or check phaseII weight

  # Get correct order of phaseI design strata vars. Must be called after tps_design_setData
  objlist$phaseI.strata <- tps_setPhaseI_design(objlist$phaseI.strata)

  # Phase I stratification
  phI_strat <- phaesI_strat.f(pcuts=objlist$pcuts, data, design.factor=objlist$phaseI.strata)

  # Phase II selection 
  # fit the working model, only adjusted for x1
  para0 <- est.H0_wei.f(data)

  phII.design <- objlist$phaseII.design.orig
  if (objlist$phaseII.design %in% c("RSD1", "RSD2", "RSD12")){
  
    # phase 2 residual-dependent sub-sampling 
    phII_sel <- try(designPhII_resid(phaseI_strat=phI_strat, n2samp=objlist$phaseII.nsamp, 
                       design=phII.design, objlist$phaseII.weight, para0$par)$sel, silent=TRUE)
    if (inherits(phII_sel, "try-error")) {
      stop("ERROR calling designPhII_resid, try changing the phaseII options")
    }
  } else {
    # phase 2 stratified sub-sampling design
    phII_sel <- design_strat.f(objlist$phaseII.nsamp, phII.design, phI_strat)
  }    

  list(phaseI=phI_strat, objlist=objlist, phaseII=phII_sel)
}

tps_setPhaseI_design <- function(strata.vars) {

  n <- length(strata.vars) 
  if (n == 1) return(strata.vars)

  ret <- NULL
  if ("x1" %in% strata.vars) {
    # Not in alphabetical order
    tlist <- list(c("x1", "d2"), c("x1", "a2", "d2"), c("d1", "d2", "x1"))
    for (i in 1:length(tlist)) {
      tmp <- tlist[[i]]
      if (all(strata.vars %in% tmp)) {
        ret <- tmp
        break
      }
    }
  } else {
    # Alphabetical order
    ret <- sort(strata.vars)
  }  
  if (!length(ret)) stop("INTERNAL CODING ERROR in tps_setPhaseI_design")
  ret
}

tps_design_setData <- function(data, objlist) {

  # Set row names
  rownames(data) <- as.character(data[, objlist$id, drop=TRUE])

  # Remove non-finite values
  vv <- c(objlist$time.lower, objlist$time.upper, 
          objlist$disease.at.baseline, objlist$disease.after.baseline,
          objlist$covariate)
  keep <- rep(TRUE, nrow(data))
  n0   <- nrow(data)
  for (v in vv) {
    tmp  <- is.finite(data[, v, drop=TRUE])
    keep <- keep & tmp
  }
  if (!all(keep)) {
    data <- data[keep, , drop=FALSE]
    nr   <- nrow(data)
    msg  <- paste0("After removing non-finite values in data, data contains ", nr, " rows")
    warning(msg)
    if (nr < 2) stop("ERROR: data contains too few rows")
  } 

  # Keep only columns we need
  vv <- c(objlist$id, objlist$time.lower, objlist$time.upper, 
          objlist$disease.at.baseline, objlist$disease.after.baseline,
          objlist$covariate)
  data <- data[, vv, drop=FALSE]

  # Rename variables
  new <- c("id", "a1", "a2", 
           "d1", "d2",
           "x1")
  colnames(data) <- new

  # Save the variable mappings
  map.origToNew         <- new
  names(map.origToNew)  <- vv
  map.newToOrig         <- vv
  names(map.newToOrig)  <- new 
  
  objlist$map.origToNew <- map.origToNew
  objlist$map.newToOrig <- map.newToOrig

  # Give objects in objlist the new names
  nms <- c("id", "time.lower", "time.upper", 
           "disease.at.baseline", "disease.after.baseline",
           "covariate", "phaseI.strata")
  for (nm in nms) {
    vars          <- objlist[[nm, exact=TRUE]]
    objlist[[nm]] <- map.origToNew[vars]
  }

  list(data=data, objlist=objlist)

}

tps_getAllPhaseII_designs <- function() {

  c("AWRSD-Smu1", "AWRSD-Smu2",
    "BRSD-eff-seq", "BRSD-eff-interac-seq",
    "FRSD-v5-eff-interac-seq",
    "RSD-Smu1", "RSD-Smu2", "RSD-Smu3",
    "RSD-Smu1-minus-2", "RSD-Smu1-plus-2", "RSD-Smu1-x1", 
    "TRSD-eff-interac-seq", "TRSD-v2-eff-interac-seq",
    "TRSD-v3-eff-interac-seq", "TRSD-v4-eff-interac-seq", 
    "TRSD-v5-eff-interac-seq",
    "WBRSD-eff-seq", "WRSD-Smu1", "WRSD-Smu2"
  )

  c("RSD1", "RSD2", "RSD12", "srs", "bal")

}
