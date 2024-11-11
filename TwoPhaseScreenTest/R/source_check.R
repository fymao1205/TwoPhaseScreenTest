check_data <- function(x, nm="data") {

  if (!is.data.frame(x)) stop(paste0("ERROR: ", nm, " must be a data frame"))
  if (nrow(x) < 2) stop(paste0("ERROR: ", nm, " contains less than 2 rows"))
  if (!ncol(x)) stop(paste0("ERROR: ", nm, " contains no columns"))
  NULL
}

check_data_est <- function(x, nm="data") {

  if (!is.data.frame(x)) stop(paste0("ERROR: ", nm, " must be a data frame"))
  if (nrow(x) < 1) stop(paste0("ERROR: ", nm, " contains no rows"))
  if (ncol(x) != 2) stop(paste0("ERROR: ", nm, " must contain two columns"))
  
  NULL
}

check_col <- function(x, nm, data, pos=0, bin=0, num=1, vals=NULL, nonNeg=0) {

  if (!isString(x)) stop(paste0("ERROR: ", nm, " must be a column name in the data"))
  valid <- colnames(data)
  if (!(x %in% valid)) {
    msg <- paste0("ERROR: ", nm, "=", getQuotedVecStr(x), " not found in the data")
    stop(msg)
  }
  vec <- data[, x, drop=TRUE]
  if (num && !is.numeric(vec)) {
    stop(paste0("ERROR: ", nm, "=", getQuotedVecStr(x), " must be numeric"))
  }

  if (pos) {
    tmp <- vec <= 0
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) stop(paste0("ERROR: ", nm, "=", getQuotedVecStr(x), " must be positive"))
  }

  if (nonNeg) {
    tmp <- vec < 0
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) stop(paste0("ERROR: ", nm, "=", getQuotedVecStr(x), " must be non-negative"))
  }

  if (bin) {
    tmp <- !is.na(vec) & !(vec %in% 0:1)
    if (any(tmp)) {
      msg <- paste0("ERROR: ", nm, "=", getQuotedVecStr(x), " must be coded 0-1")
      stop(msg)
    }
  }

  if (length(vals)) {
    tmp <- is.finite(vec) & !(vec %in% vals)
    if (any(tmp)) {
      str <- paste0(vals, collapse=", ") 
      msg <- paste0("ERROR: ", nm, "=", getQuotedVecStr(x), " must be coded ", str)
      stop(msg)
    }
  }


  NULL
}

check_cols <- function(x, nm, data, n=0, min=0, num=0) {

  len <- length(x)
  if (min && (len < min)) {
    msg <- paste0("ERROR: ", nm, " must be a vector of length >= ", min)
    stop(msg) 
  }
  if (len) {
    if (n && (n != len)) {
      msg <- paste0("ERROR: ", nm, " must be a vector of length ", n)
      stop(msg) 
    }
    if (!is.vector(x) || !is.character(x)) {
      msg <- paste0("ERROR: ", nm, " must be a vector of column names in the data")
      stop(msg) 
    }
    for (v in x) check_col(v, nm, data, pos=0, bin=0, num=num) 
    if (len != length(unique(x))) {
      msg <- paste0("ERROR: ", nm, " does not contain unique column names")
      stop(msg) 
    }
  }
  NULL
}

check_strata.cols <- function(x, data, valid, nm="phase1.strata") {

  if (!length(x)) x <- valid[1]
  x <- unique(x)
  check_cols(x, nm, data, n=0, min=0, num=1)
  tmp <- !(x %in% valid) 
  if (any(tmp)) {
    str <- getQuotedVecStr(x[tmp])
    msg <- paste0("ERROR: ", nm, " contains the invalid column(s) ", str)
    stop(msg)
  } 
  x
}

check_time.col <- function(x, data, nm="time.col") {

  check_col(x, nm, data, pos=1, bin=0)
  NULL
}

check_status.col <- function(x, data, nm="status.col") {

  check_col(x, nm, data, pos=0, bin=1)
  NULL
}

check_exposure.col <- function(x, data, nm="exposure.col") {

  check_col(x, nm, data, pos=0, bin=0)
  NULL
}



getQuotedVecStr <- function(x, sep=",") {
  ret <- paste0("'", x, "'")
  ret <- paste0(ret, collapse=sep)
  ret
}

isString <- function(x) {
  (length(x) == 1) && is.character(x)
}

check_num <- function(x, name, len=1, pos=1, min=NULL, max=NULL) {

  xlen <- length(x)
  if (!len && !xlen) return(NULL)
  if (len && (xlen != len)) {
    stop(paste0("ERROR: ", name, " must be a numeric value"))
  }
  tmp <- !is.finite(x)
  if (any(tmp)) stop(paste0("ERROR: ", name, " must be a numeric value"))
  if (pos && (any(x <= 0))) {
    stop(paste0("ERROR: ", name, " must be a positive value"))
  }

  if (length(min) && any(x < min)) {
    stop(paste0("ERROR: ", name, " must be >= ", min))
  }
  if (length(max) && any(x > max)) {
    stop(paste0("ERROR: ", name, " must be <= ", max))
  }

  NULL
}

check_int <- function(x, nm, len=1, pos=1, min=NULL, max=NULL) {

  check_num(x, nm, len=len, pos=pos, min=min, max=max)
  if (x != floor(x)) {
    msg <- paste0("ERROR: ", nm, " must be an integer")
    stop(msg)
  }
  NULL
}

check_log <- function(x, name) {

  if (length(x) != 1) stop(paste0("ERROR: ", name, " must be TRUE or FALSE"))
  tmp <- x %in% c(TRUE, FALSE)
  if (!tmp) stop(paste0("ERROR: ", name, " must be TRUE or FALSE"))
  NULL
}

check_str <- function(x, nm, valid) {
  
  if (!isString(x)) {
    err <- getQuotedVecStr(valid, sep=",")
    msg <- paste0("ERROR: ", nm, " must be one of ", err)
    stop(msg)
  }
  if (!(x %in% valid)) {
    err <- getQuotedVecStr(valid, sep=",")
    msg <- paste0("ERROR: ", nm, " must be one of ", err)
    stop(msg)
  }
  NULL
}

isList <- function(x) {

  ret <- is.list(x) && ("list" %in% class(x))
  ret
}

checkRequiredListNames <- function(x, req, name) {

  tmp  <- !(req %in% names(x))
  miss <- req[tmp]
  if (length(miss)) {
    tmp <- paste0(getQuotedVecStr(miss) , collapse=", ")
    msg <- paste0("ERROR: the objects ", tmp, " not found in ", name)
    stop(msg)  
  }
  
  NULL

} # END: checkRequiredListNames

check.list <- function(x, name, req) {

  if (!isList(x)) stop(paste0("ERROR: ", name, " must be a list"))
  if (!length(x)) stop(paste0("ERROR: ", name, " is empty"))
  checkRequiredListNames(x, req, name)

  NULL 

} # END: check.list

check_phaseII.design <- function(x) {

  valid <- tps_getAllPhaseII_designs()
  check_str(x, "phaseII.design", valid)
  NULL
}

check_n2.samp <- function(x, data, nm="n2.samp") {

  if (length(x)) {
    check_int(x, nm, len=1, pos=1, min=1, max=nrow(data)-1)
  } 
  x
}

check_scenario <- function(x) {

  valid <- c("a", "b", "c")
  check_str(x, "scenario", valid)
  NULL
}

check_exposure.type <- function(x) {

  nm    <- "exposure.type"
  if (length(x)) {
    valid <- c("binary", "cont")
    check_str(x, nm, valid)
  }
  NULL
}

check_cuts <- function(x) {

  nm  <- "pcuts"
  def <- c(1/3, 2/3)
  if (!length(x)) x <- def
  if (!is.numeric(x) || !is.vector(x)) {
    stop(paste0("ERROR: ", nm, " must be a numeric vector"))
  }
  tmp <- (x < 0) | (x > 1) | !is.finite(x)
  tmp[is.na(tmp)] <- TRUE
  if (any(tmp)) {
    stop(paste0("ERROR: ", nm, " must have values in [0, 1]"))
  }

  x
}

check_design.obj <- function(x, nm="design.obj") {

  req <- c("selected.ids", "phaseI.strata.freq", "phaseII.strata.freq",
           "data", "options", "variable.map")
  check.list(x, nm, req) 

  NULL
}

check_est.obj <- function(x, nm="est.obj") {

  req <- c("est", "cov", "brks", "interac.ind")
  check.list(x, nm, req) 

  NULL
}


check_interac.ind <- function(x, nm="interac.ind") {

  if (length(x) != 2) stop(paste0("ERROR: ", nm, " must be a numeric vector of length 2"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a numeric vector of length 2"))
  if (!is.vector(x))  stop(paste0("ERROR: ", nm, " must be a numeric vector of length 2"))
  tmp <- all(x %in% 0:1)
  if (!tmp) {
    msg <- paste0("ERROR: ", nm, " must be one of c(0,0), c(0,1), c(1,0) or c(1,1)")
    stop(msg)
  }
  NULL
}
