
# different quantities measuring the diagnostic accuracies ...


cumrisk_g_x2_pwc.f <- function(t, x2, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  #res1 <- exp(logL_cal_pwc(2, 1, 0, t, 1, x2, gam, alp, beta, eta, p1, brks))#*p1
  #res0 <- exp(logL_cal_pwc(2, 1, 0, t, 0, x2, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res1 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 1, x2, gam, alp, beta, eta, p1, brks))#*p1
  res0 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 0, x2, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  p11 = expit.f(sum(eta))*p1
  p10 = expit.f(eta[1])*(1-p1)
  
  denom=x2*(p11+p10)+(1-x2)*( 1-p11-p10 )
  
  return((res0+res1)/denom)
}

ase.cumrisk_g_x2_pwc.f <- function(grad, para.varmat, t, x2, para, brks, interac.ind){
  
  res0 <- cumrisk_g_x2_pwc.f(t, x2, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- cumrisk_g_x2_pwc.f(t, x2, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, partial=dres, ase=ase, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

cumrisk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  #res1 <- exp(logL_cal_pwc(2, 1, 0, t, x1, 1, gam, alp, beta, eta, p1, brks))#*p1
  #res0 <- exp(logL_cal_pwc(2, 1, 0, t, x1, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res1 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, x1, 1, gam, alp, beta, eta, p1, brks))#*p1
  res0 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, x1, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)

  
  denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1)/denom)
}

ase.cumrisk_g_x1_pwc.f <- function(grad, para.varmat, t, x1, para, brks, interac.ind){
  
  res0 <- cumrisk_g_x1_pwc.f(t, x1, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- cumrisk_g_x1_pwc.f(t, x1, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, ase=ase, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

pre1_risk_g_x2_pwc.f <- function(x2, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  #alp = exp(para[2+1:len])
  #beta = c(para[2+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x2_1 = (x2*p11+(1-x2)*(p1-p11))
  p_x2_0 = (x2*p10+(1-x2)*(1-p1-p10))
  
  res1 <- p_x2_1*expit.f(gam[1]+gam[2]+gam[3]*x2)#exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- p_x2_0*expit.f(gam[1]+gam[3]*x2)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x2*(p11+p10)+(1-x2)*( 1-p11-p10 )
  
  return((res0+res1)/denom)
}

pre1_risk_g_x1_pwc.f <- function(x1, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  #alp = exp(para[2+1:len])
  #beta = c(para[2+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  
  
  p1_x1=expit.f(eta[1]+eta[2]*x1)#*p1
  #p_x1_0=expit.f(eta[1])#*(1-p1)
  
  #p_x1_1 = (x1*p11+(1-x1)*p10)
  #p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  res1 <- p1_x1*expit.f(gam[1]+gam[2]*x1+gam[3]+gam[4]*x1)#exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- (1-p1_x1)*expit.f(gam[1]+gam[2]*x1)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  #denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1))
}

ase.pre1_risk_g_x1_pwc.f <- function(grad, para.varmat, x1, para, brks, interac.ind){
  
  res0 <- pre1_risk_g_x1_pwc.f(x1, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- pre1_risk_g_x1_pwc.f(x1, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

ase.pre1_risk_g_x2_pwc.f <- function(grad, para.varmat, x2, para, brks, interac.ind){
  
  res0 <- pre1_risk_g_x2_pwc.f(x2, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- pre1_risk_g_x2_pwc.f(x2, para.p, brks, interac.ind)
    (res0-res0.p)/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=res0, lb=res0-1.96*ase, ub=res0+1.96*ase))
}

inc1_risk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  #p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  #gam = c(para[2+len+len.beta+1:len.gam],0)
  
  p1_x1=expit.f(eta[1]+eta[2]*x1)#*p1
  #p10=expit.f(eta[1])*(1-p1)
  
  #p_x1_1 = (x1*p11+(1-x1)*p10)
  #p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  levs_x1_1 = alp*exp(beta[1]*x1+beta[2]+beta[3]*x1); 
  levs_x1_0 = alp*exp(beta[1]*x1); 
  surv_x1_1 = RtoC_ppwc(t, brks, levs_x1_1, lower=FALSE,0)
  surv_x1_0 = RtoC_ppwc(t, brks, levs_x1_0, lower=FALSE,0)
  
  res1 <- p1_x1 * (1-surv_x1_1) #exp(logL_cal_pwc(0, NA, 1, t, x1, 1, gam, alp, beta, eta, brks))*p1
  res0 <- (1-p1_x1) * (1-surv_x1_0) #exp(logL_cal_pwc(0, NA, 1, t, x1, 0, gam, alp, beta, eta, brks))*(1-p1)
  
  #denom=x1*p1+(1-x1)*( 1-p1)
  
  #return((res0+res1)/denom)
  return(res0+res1)
}

inc1_risk_g_x2_pwc.f <- function(t, x2, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  #gam = c(para[2+len+len.beta+1:len.gam],0)
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x2_1 = (x2*p11+(1-x2)*(p1-p11))
  p_x2_0 = (x2*p10+(1-x2)*(1-p1-p10))
  
  levs_x2_1 = alp*exp(beta[1]+beta[2]*x2+beta[3]*x2); 
  levs_x2_0 = alp*exp(beta[2]*x2); 
  surv_x2_1 = RtoC_ppwc(t, brks, levs_x2_1, lower=FALSE,0)
  surv_x2_0 = RtoC_ppwc(t, brks, levs_x2_0, lower=FALSE,0)
  
  res1 <- p_x2_1 * (1-surv_x2_1) #exp(logL_cal_pwc(0, NA, 1, t, x1, 1, gam, alp, beta, eta, brks))*p1
  res0 <- p_x2_0 * (1-surv_x2_0) #exp(logL_cal_pwc(0, NA, 1, t, x1, 0, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x2*(p11+p10)+(1-x2)*( 1-p11-p10 )
  
  return((res0+res1)/denom)
}

wgt1_risk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x1_1 = (x1*p11+(1-x1)*p10)
  p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  
  levs_x1_1 = alp*exp(beta[1]*x1+beta[2]+beta[3]*x1); 
  levs_x1_0 = alp*exp(beta[1]*x1); 
  surv_x1_1 = RtoC_ppwc(t, brks, levs_x1_1, lower=FALSE,0)
  surv_x1_0 = RtoC_ppwc(t, brks, levs_x1_0, lower=FALSE,0)
  
  pz_x1_1 = expit.f(gam[1]+gam[2]*x1+gam[3]+gam[4]*x1)
  pz_x1_0 = expit.f(gam[1]+gam[2]*x1)
  
  res1 <- p_x1_1*pz_x1_1/(1-(1-pz_x1_1)*surv_x1_1)  #exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- p_x1_0*pz_x1_0/(1-(1-pz_x1_0)*surv_x1_0)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1)/denom)
  
  #return(exp(res0))
  
}

wgt2_risk_g_x1_pwc.f <- function(t, x1, para, brks, interac.ind){
  
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  p11=expit.f(sum(eta))*p1
  p10=expit.f(eta[1])*(1-p1)
  
  p_x1_1 = (x1*p11+(1-x1)*p10)
  p_x1_0 = (x1*(p1-p11)+(1-x1)*(1-p1-p10))
  
  
  levs_x1_1 = alp*exp(beta[1]*x1+beta[2]+beta[3]*x1); 
  levs_x1_0 = alp*exp(beta[1]*x1); 
  surv_x1_1 = RtoC_ppwc(t, brks, levs_x1_1, lower=FALSE,0)
  surv_x1_0 = RtoC_ppwc(t, brks, levs_x1_0, lower=FALSE,0)
  
  pz_x1_1 = expit.f(gam[1]+gam[2]*x1+gam[3]+gam[4]*x1)
  pz_x1_0 = expit.f(gam[1]+gam[2]*x1)
  
  res1 <- p_x1_1*(1-pz_x1_1)*(1-surv_x1_1)/(1-(1-pz_x1_1)*surv_x1_1)  #exp(logL_cal_pwc(0, NA, 1, t, 1, x2, gam, alp, beta, eta, brks))*p1
  res0 <- p_x1_0*(1-pz_x1_0)*(1-surv_x1_0)/(1-(1-pz_x1_0)*surv_x1_0)#exp(logL_cal_pwc(0, NA, 1, t, 0, x2, gam, alp, beta, eta, brks))*(1-p1)
  
  denom=x1*p1+(1-x1)*( 1-p1)
  
  return((res0+res1)/denom)
  
  #return(exp(res0))
  
}


# TPR=P(X1=1|D0=1); FPR=P(X1=1|D0=0)
tpr.fpr.auc_of_x1.f <- function(para, brks, interac.ind){
  
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  
  # P(x1, x2)
  p11 = expit.f( sum(eta) )*p1
  p10 = p1 - p11
  p01 = expit.f( eta[1] )*(1-p1)
  p00 = (1-p1)-p01
  
  # numerator
  res11 = expit.f( cbind(1, 1, 1, 1) %*% gam )*p11 # P(D0=1|x=(1,1))*P(x=(1,1))
  res10 = expit.f( cbind(1, 1, 0, 0) %*% gam )*p10 # P(D0=1|x=(1,0))*P(x=(1,0))
  res01 = expit.f( cbind(1, 0, 1, 0) %*% gam )*p01 # P(D0=1|x=(0,1))*P(x=(0,1))
  res00 = expit.f( cbind(1, 0, 0, 0) %*% gam )*p00 # P(D0=1|x=(0,0))*P(x=(0,0))
  denom1 = res11 + res10 + res01 + res00 # P(D0=1)
  
  tpr= (res11+res10)/denom1
  fpr= (p11+p10-(res11+res10))/(1-denom1)
  
  
  auc=tpr*(1-fpr)
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

# TPR=P(X2=1|D0=1); FPR=P(X2=1|D0=0)
tpr.fpr.auc_of_x2.f <- function(para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  
  # P(x1, x2)
  p11 = expit.f( sum(eta) )*p1
  p10 = p1 - p11
  p01 = expit.f( eta[1] )*(1-p1)
  p00 = (1-p1)-p01
  
  # numerator
  res11 = expit.f( cbind(1, 1, 1, 1) %*% gam )*p11 # P(D0=1|x=(1,1))*P(x=(1,1))
  res10 = expit.f( cbind(1, 1, 0, 0) %*% gam )*p10 # P(D0=1|x=(1,0))*P(x=(1,0))
  res01 = expit.f( cbind(1, 0, 1, 0) %*% gam )*p01 # P(D0=1|x=(0,1))*P(x=(0,1))
  res00 = expit.f( cbind(1, 0, 0, 0) %*% gam )*p00 # P(D0=1|x=(0,0))*P(x=(0,0))
  denom1 = res11 + res10 + res01 + res00 # P(D0=1)
  
  tpr= (res11+res01)/denom1
  fpr= (p11+p01-(res11+res01))/(1-denom1)
  
  auc=tpr*(1-fpr)
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

# TPR=P(X2=1|D(t)=1); FPR=P(X2=1|D(t)=0)
dyn.tpr.fpr.auc_of_x2.f <- function(t, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  #res1 <- exp(logL_cal_pwc(2, 1, 0, t, 1, x2, gam, alp, beta, eta, p1, brks))#*p1
  #res0 <- exp(logL_cal_pwc(2, 1, 0, t, 0, x2, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  #res11 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 1, gam, alp, beta, eta, p1, brks))#*p1
  #res01 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 1, gam, alp, beta, eta, p1, brks))#*(1-p1)
  #res10 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 0, gam, alp, beta, eta, p1, brks))#*p1
  #res00 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res11 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 1, 1, gam, alp, beta, eta, p1, brks))#*p1
  res01 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 0, 1, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res10 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 1, 0, gam, alp, beta, eta, p1, brks))#*p1
  res00 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 0, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  denom1 = res11 + res01 + res10 + res00
  
  p11 = expit.f(sum(eta))*p1
  p10 = expit.f(eta[1])*(1-p1)
  
  tpr= (res11+res01)/denom1
  fpr= (p11+p10-(res11+res01))/(1-denom1)
  
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

# TPR=P(X2=1|D(t)=1); FPR=P(X2=1|D(t)=0)
dyn.tpr.fpr.auc_of_x1.f <- function(t, para, brks, interac.ind){
  
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
  
  eta = para[1:2]
  p1 = expit.f(para[3])
  alp = exp(para[3+1:len])
  beta = c(para[3+len+1:len.beta], 0)
  gam = c(para[3+len+len.beta+1:len.gam],0)
  #p1 = para[2+len+len.beta+len.gam+1]
  
  #res1 <- exp(logL_cal_pwc(2, 1, 0, t, 1, x2, gam, alp, beta, eta, p1, brks))#*p1
  #res0 <- exp(logL_cal_pwc(2, 1, 0, t, 0, x2, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  #res11 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 1, gam, alp, beta, eta, p1, brks))#*p1
  #res01 <- exp(logL_cal_pwc(2, 1, 0, t, 1, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  #res10 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 1, gam, alp, beta, eta, p1, brks))#*p1
  #res00 <- exp(logL_cal_pwc(2, 1, 0, t, 0, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res11 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 1, 1, gam, alp, beta, eta, p1, brks))#*p1
  res01 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 1, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  res10 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 0, 1, gam, alp, beta, eta, p1, brks))#*p1
  res00 <- exp(RtoC_logL_cal_pwc(2, 1, 0, t, 0, 0, gam, alp, beta, eta, p1, brks))#*(1-p1)
  
  denom1 = res11 + res01 + res10 + res00
  
  tpr= (res11+res01)/denom1
  fpr= (p1-(res11+res01))/(1-denom1)
  
  auc.ties=0.5*(tpr+(1-fpr))
  
  return(list(tpr=tpr, fpr=fpr, #auc=auc, 
              auc.ties=auc.ties))
}

ase.tpr.fpr.auc_of_x2.f <- function(grad, para.varmat, para, brks, interac.ind){
  
  res0 <- tpr.fpr.auc_of_x2.f(para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- tpr.fpr.auc_of_x2.f(para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.tpr.fpr.auc_of_x1.f <- function(grad, para.varmat, para, brks, interac.ind){
  
  res0 <- tpr.fpr.auc_of_x1.f(para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- tpr.fpr.auc_of_x1.f(para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.dyn.tpr.fpr.auc_of_x2.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x2.f(t, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x2.f(t, para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.dyn.tpr.fpr.auc_of_x1.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x1.f(t, para, brks, interac.ind)
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x1.f(t, para.p, brks, interac.ind)
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  
  
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}


ase.dyn.auc_of_x2.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x2.f(t, para, brks, interac.ind)$auc.ties
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x2.f(t, para.p, brks, interac.ind)$auc.ties
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}

ase.dyn.auc_of_x1.f <- function(t, grad, para.varmat, para, brks, interac.ind){
  
  res0 <- dyn.tpr.fpr.auc_of_x1.f(t, para, brks, interac.ind)$auc.ties
  
  dres.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    res0.p <- dyn.tpr.fpr.auc_of_x1.f(t, para.p, brks, interac.ind)$auc.ties
    (as.vector(unlist(res0))-as.vector(unlist(res0.p)))/grad 
  })
  
  dres <- do.call("cbind", dres.list)
  
  vcov.delta=dres %*% para.varmat %*% t(dres)
  
  ase=sqrt(diag(vcov.delta))
  return(list(est=as.vector(unlist(res0)), 
              ase=ase,
              lb=as.vector(unlist(res0))-1.96*ase, 
              ub=as.vector(unlist(res0))+1.96*ase))
}





