
wgt.d1_g_obs.f <- function(z, d1, d2, t1, t2, x1, gam, alp, beta){

  piZ = expit.f(gam[1]+gam[2]*x1)
  expcovT=exp(beta[1]*x1)
  Hft1 <- (alp[1]*t1)^(alp[2])*expcovT 
  survft1 <- exp(-Hft1)
  Hft2 <- (alp[1]*t2)^(alp[2])*expcovT 
  survft2 <- exp(-Hft2)
  
  res1 <- ifelse(d2==0, 0, ifelse(d1>1, piZ/( 1- (1-piZ)*survft2 ), d1))
  
  res <- z*res1+(1-z)*(1-res1) 

  return(res)
}

Smu_vec.f <- function(d1, d2, t1, t2, x1, para.H0){
  
  alp=c(exp(para.H0[1]), exp(para.H0[2])); beta=para.H0[3]; gam=para.H0[4:5]
  piZ = expit.f(gam[1]+gam[2]*x1)
  expcovT=exp(beta*x1)
  Hft1 <- (alp[1]*t1)^(alp[2])*expcovT #exp(-exp(alp+beta*v)*t) #vppwc(t, brks, exp(alp), 0, 0)
  survft1 <- exp(-Hft1)
  Hft2 <- (alp[1]*t2)^(alp[2])*expcovT #exp(-exp(alp+beta*v)*t) #vppwc(t, brks, exp(alp), 0, 0)
  survft2 <- exp(-Hft2)
  
  wgt <- wgt.d1_g_obs.f(1, d1, d2, t1, t2, x1, gam, alp, beta)

  res1 <- wgt-piZ 
  
  res20 <- ifelse(d1==1, 0, survft2*(-Hft2))
  res2 <- ifelse( d1>1,  
                  ifelse(d2==0, 0, (1-piZ)*(1-survft2)/(1-(1-piZ)*survft2)*Hft2),
                  (1-d1)*( d2*(  (t1>0)*survft1*(-Hft1)) - res20 )/(d2*survft1-survft2)
                  )
  
  res <- list(Smu1=res1, Smu2=res2, Smu3=res1-res2)
  
  return(res)
  
}

logLn_cal1_mu_wei.f <- function(mu1, mu2, d1, d2, t1, t2, alp){

  expcovZ= exp(mu1)
  expcovT=exp(mu2)
  Hft1 = (alp[1]*t1)^(alp[2])*expcovT
  Sft1 = exp(-Hft1)
  Hft2 = (alp[1]*t2)^(alp[2])*expcovT
  Sft2 = exp(-Hft2)
  piZ = expcovZ/(1+expcovZ)
  
  res <- log(piZ)
  res[d1>1 & d2==1] <- (log((1-Sft2)*(1-piZ)+piZ))[d1>1 & d2==1]
  res[d1>1 & d2==0] <- rep(0, sum(d1>1 & d2==0))
  res[d2>1] <- rep(0, sum(d2>1))
  res[d1==0 & d2==1] <- (log( Sft1-Sft2 )+log(1-piZ))[d1==0 & d2==1]
  res[d1==0 & d2==0] <- (log(1-piZ)-Hft2)[d1==0 & d2==0]
    
  return(res)
}


Smu_num_vec.f <- function(grad, d1, d2, a1, a2, x1, para.H0){
  
  alp=c(exp(para.H0[1]), exp(para.H0[2])); beta=para.H0[3]; gam=para.H0[4:5]
  mu1=gam[1]+gam[2]*x1
  mu2=beta[1]*x1
  
  l0 = logLn_cal1_mu_wei.f(mu1, mu2, d1, d2, a1, a2, alp)
  
  lp_mu2 = logLn_cal1_mu_wei.f(mu1, mu2+grad, d1, d2, a1, a2, alp)
  
  lp_mu1 = logLn_cal1_mu_wei.f(mu1+grad, mu2, d1, d2, a1, a2, alp)
  
  res <- list(Smu1=(lp_mu1-l0)/grad, Smu2=(lp_mu2-l0)/grad)
  
  return(res)
}

est.H0_wei.f <- function(dat){
  
  obj.f <- function(para){
    
    alp = exp(para[1:2])
    beta = c(para[3],0,0)
    gamma = c(para[4:5],0,0)
  
    #res0 <- log.calLn.f(dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, alp, beta, gamma)
    #res0 <- logLn_cal1_wei(dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gamma, alp, beta)
    res0 <- RtoC_logLn_cal1_wei(dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gamma, alp, beta)
 
    res <- -sum(res0)
    
    # print(para)
    # print(res)
    
    return(res)
    
  }
  
  est <- optim(par=rep(0.01, 5), obj.f, method=("L-BFGS-B"), hessian=TRUE)
  #est <- nlm(f=obj.f, p=rep(0.01, p1+p2+p3+4), hessian=TRUE)
  return(est)
}

designPhII_resid <- function(phaseI_strat, n2samp, design="RSD-Smu1", w=NULL,
                                 para.H0){
  
  # tp prevent compiler warning
  R <- NULL

  # qzi <- Smu_vec.f(phaseI_strat$dt_ext$d1, phaseI_strat$dt_ext$d2,
  #                  phaseI_strat$dt_ext$a1, phaseI_strat$dt_ext$a2,
  #                  phaseI_strat$dt_ext$x1, para.H0)
  
  qzi <- Smu_num_vec.f(1e-06, phaseI_strat$dt_ext$d1, phaseI_strat$dt_ext$d2,
                        phaseI_strat$dt_ext$a1, phaseI_strat$dt_ext$a2,
                        phaseI_strat$dt_ext$x1, para.H0)
  
  k= round(min(n2samp/2, sum(phaseI_strat$dt_ext_b$x1==0)/2, sum(phaseI_strat$dt_ext_b$x1==1)/2)*0.15)

  # Set min of k to 10 (email 2024-10-25)
  if (k < 10) k <- 10

  #k=round(n2samp/2*0.05)
  if(design %in% c("RSD-Smu1")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[1]], phaseI_strat) 
  }else if(design %in% c("RSD-Smu1-minus-2")){
    rsd=qzi[[1]]-qzi[[2]]
    res0 <- rkd_ext_opt.f(k, n2samp, rsd, phaseI_strat)
  }else if(design %in% c("RSD-Smu1-plus-2")){
    rsd=qzi[[1]]+qzi[[2]]
    res0 <- rkd_ext_opt.f(k, n2samp, rsd, phaseI_strat)
  }else if(design %in% c("RSD-Smu1-x1")){
    rsd=qzi[[1]]*phaseI_strat$dt_ext$x1
    res0 <- rkd_ext_opt.f(k, n2samp, rsd, phaseI_strat)
  }else if(design %in% c("RSD-Smu2")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[2]], phaseI_strat)
  }else if(design %in% c("RSD-Smu3")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[3]], phaseI_strat)
  }else if(design %in% c("WRSD-Smu1")){
    
    # covariate model
    fit_x2_g_x1 <- glm(x2~x1, family=binomial (link=logit), data=phaseI_strat$dt_ext)
    sd_x2_g_x1 <- sqrt(fit_x2_g_x1$fitted.values*(1-fit_x2_g_x1$fitted.values))
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[1]]*sd_x2_g_x1, phaseI_strat)
    
  }else if(design %in% c("WRSD-Smu2")){
    
    # covariate model
    fit_x2_g_x1 <- glm(x2~x1, family=binomial (link=logit), data=phaseI_strat$dt_ext)
    sd_x2_g_x1 <- sqrt(fit_x2_g_x1$fitted.values*(1-fit_x2_g_x1$fitted.values))
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[2]]*sd_x2_g_x1, phaseI_strat)
    
  }else if(design %in% "AWRSD-Smu1"){

      n2a=round(n2samp*w)
      ka=round(n2a/2*0.15)
      #res01 <- rkd_ext_opt.f(ka, n2a, Seff1, phaseI_strat)
      res01 <- design_strat.f(n2a, design="srs", phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
      
      # covariate model
      fit_x2_g_x1 <- glm(x2~x1, family=binomial (link=logit), data=subset(dt_a_use, R==1))
      p1_x2_g_x1 <- predict(fit_x2_g_x1, newdata = dt_a_use, type="response")
      sd_x2_g_x1 <- sqrt(p1_x2_g_x1*(1-p1_x2_g_x1))
      
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, qzi[[1]]*sd_x2_g_x1, phaseI_strat, res01$s_id)
    
  }else if(design %in% "AWRSD-Smu2"){
    
    n2a=round(n2samp*w)
    #ka=round(n2a/2*0.15)
    ka= round(min(n2a/2, sum(phaseI_strat$dt_ext$x1==0)/2, sum(phaseI_strat$dt_ext$x1==1)/2)*0.15)
    
    #res01 <- rkd_ext_opt.f(ka, n2a, Seff1, phaseI_strat)
    res01 <- design_strat.f(n2a, design="srs", phaseI_strat)
    
    dt_a_use <- phaseI_strat$dt_ext #dt_c
    dt_a_use$R<-1
    dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
    
    # covariate model
    fit_x2_g_x1 <- glm(x2~x1, family=binomial (link=logit), data=subset(dt_a_use, R==1))
    
    dt_b=subset(dt_a_use, R==0)
    
    p1_x2_g_x1 <- predict(fit_x2_g_x1, newdata = dt_a_use, type="response")
    sd_x2_g_x1 <- sqrt(p1_x2_g_x1*(1-p1_x2_g_x1))
    
    n2b=n2samp-n2a
    #kb=round(n2b/2*0.15)
    kb= round(min(n2b/2, sum(dt_b$x1==0)/2, sum(dt_b$x1==1)/2)*0.15)
    phaseI_strat.b=phaseI_strat
    phaseI_strat.b$dt_ext=dt_b
    phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
    res0 <- rkd_ext_opt_2st.f(kb, n2b, qzi[[2]]*sd_x2_g_x1, phaseI_strat, res01$s_id)
    
  }else if(design %in% "BRSD-eff-seq"){
    
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    covS=cov(mat)
    Seff1=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[1,2]/covS[1,1]*mat[,1]
    
    
    if(w==0 | w==1){
      qzi_use=w*Seff1+(1-w)*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a=round(n2samp*w)
      #ka=round(n2a/2*0.15)
      ka= round(min(n2a/2, sum(phaseI_strat$dt_ext$x1==0)/2, sum(phaseI_strat$dt_ext$x1==1)/2)*0.15)
      
      # Set min value to 10 (email 2024-10-25). 
      if (ka < 10) ka <- 10
      
      res01 <- rkd_ext_opt.f(ka, n2a, Seff1, phaseI_strat)
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a
      #kb=round(n2b/2*0.15)
      kb= round(min(n2b/2, sum(dt_b$x1==0)/2, sum(dt_b$x1==1)/2)*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      
    }
    
  }else if(design %in% "WBRSD-eff-seq"){
    
    # covariate model
    fit_x2_g_x1 <- glm(x2~x1, family=binomial (link=logit), data=phaseI_strat$dt_ext)
    sd_x2_g_x1 <- sqrt(fit_x2_g_x1$fitted.values*(1-fit_x2_g_x1$fitted.values))
    
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    covS=cov(mat)
    Seff1=(mat[,1]-covS[2,1]/covS[2,2]*mat[,2])*sd_x2_g_x1
    Seff2=(mat[,2]-covS[1,2]/covS[1,1]*mat[,1])*sd_x2_g_x1
    
    if(w==0 | w==1){
      qzi_use=w*Seff1+(1-w)*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a=round(n2samp*w)
      ka=round(n2a/2*0.15)
      res01 <- rkd_ext_opt.f(ka, n2a, Seff1, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      
    }
    
  }else if(design %in% "BRSD-eff-interac-seq"){
    mat=as.matrix(cbind(qzi[[1]]*(1+phaseI_strat$dt_ext$x1), qzi[[2]]))
    
    covS=cov(mat)
    #Seff1.main=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    #Seff1.interac=mat[,2]--covS[1,2]/covS[1,1]*mat[,1]-covS[3,1]/covS[3,3]*mat[,3]
    
    Seff1=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[1,2]/covS[1,1]*mat[,1]
    
    if(w==0 | w==1){
      qzi_use=w*Seff1+(1-w)*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a=round(n2samp*w)
      ka=round(n2a/2*0.15)
      res01 <- rkd_ext_opt.f(ka, n2a, Seff1, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      
    }
    
  }else if(design %in% "TRSD-eff-interac-seq"){
    mat=as.matrix(cbind(qzi[[1]], qzi[[1]]*phaseI_strat$dt_ext$x1, qzi[[2]]))
    
    covS=cov(mat)
    Seff1.main=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    Seff1.interac=mat[,2]-covS[1,2]/covS[1,1]*mat[,1]-covS[3,1]/covS[3,3]*mat[,3]
    
    #Seff1=mat[,1]+mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    
    Seff2=mat[,3]-covS[1,3]/covS[1,1]*mat[,1]-covS[2,3]/covS[2,2]*mat[,2]
    
    if(sum(w)==0 | sum(w)==1){
      qzi_use=w[1]*Seff1.main+w[2]*Seff1.interac+(1-sum(w))*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a.main=round(n2samp*w[1])
      ka.main=round(n2a.main/2*0.15)
      res01.main <- rkd_ext_opt.f(ka.main, n2a.main, Seff1.main, phaseI_strat)
      
      n2a.interac=round(n2samp*w[2])
      ka.interac=round(n2a.interac/2*0.15)
      res01.interac <- rkd_ext_opt.f(ka.interac, n2a.interac, Seff1.interac, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% c(res01.main$s_id, res01.interac$s_id)), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a.main-n2a.interac
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01.main$s_m-res01.interac$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, c(res01.main$s_id, res01.interac$s_id))
      
    }
    
  }else if(design %in% "TRSD-v2-eff-interac-seq"){
    mat=as.matrix(cbind(qzi[[1]], qzi[[1]]*phaseI_strat$dt_ext$x1, qzi[[2]]))
    
    covS=cov(mat)
    Seff1.main=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    Seff1.interac=mat[,2] #-covS[3,2]/covS[3,3]*mat[,3]
    
    #Seff1=mat[,1]+mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    
    Seff2=mat[,3]-covS[1,3]/covS[1,1]*mat[,1]-covS[2,3]/covS[2,2]*mat[,2]
    
    if(sum(w)==0 | sum(w)==1){
      qzi_use=w[1]*Seff1.main+w[2]*Seff1.interac+(1-sum(w))*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a.main=round(n2samp*w[1])
      ka.main=round(n2a.main/2*0.15)
      res01.main <- rkd_ext_opt.f(ka.main, n2a.main, Seff1.main, phaseI_strat)
      
      n2a.interac=round(n2samp*w[2])
      ka.interac=round(n2a.interac/2*0.15)
      res01.interac <- rkd_ext_opt.f(ka.interac, n2a.interac, Seff1.interac, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% c(res01.main$s_id, res01.interac$s_id)), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a.main-n2a.interac
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01.main$s_m-res01.interac$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, c(res01.main$s_id, res01.interac$s_id))
      
    }
    
  }else if(design %in% "TRSD-v3-eff-interac-seq"){
    mat=as.matrix(cbind(qzi[[1]], qzi[[1]]*phaseI_strat$dt_ext$x1, qzi[[2]]))
    
    covS=cov(mat)
    Seff1.main=mat[,1]-covS[3,1]/covS[3,3]*mat[,3]
    #Seff1.interac=mat[,2]-covS[3,2]/covS[3,3]*mat[,3]
    Seff1.interac=mat[,2]-covS[3,1]/covS[3,3]*mat[,3] #-covS[1,2]/covS[1,1]*mat[,1]
    
    #Seff1=mat[,1]+mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    
    Seff2=mat[,3]-covS[1,3]/covS[1,1]*mat[,1] #-covS[2,3]/covS[2,2]*mat[,2]
    
    if(sum(w)==0 | sum(w)==1){
      qzi_use=w[1]*Seff1.main+w[2]*Seff1.interac+(1-sum(w))*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a.main=round(n2samp*w[1])
      ka.main=round(n2a.main/2*0.15)
      res01.main <- rkd_ext_opt.f(ka.main, n2a.main, Seff1.main, phaseI_strat)
      
      n2a.interac=round(n2samp*w[2])
      ka.interac=round(n2a.interac/2*0.15)
      res01.interac <- rkd_ext_opt.f(ka.interac, n2a.interac, Seff1.interac, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% c(res01.main$s_id, res01.interac$s_id)), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a.main-n2a.interac
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01.main$s_m-res01.interac$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, c(res01.main$s_id, res01.interac$s_id))
      
    }
    
  }else if(design %in% "TRSD-v4-eff-interac-seq"){
    mat=as.matrix(cbind(qzi[[1]], qzi[[1]]*phaseI_strat$dt_ext$x1, qzi[[2]]))
    
    covS=cov(mat)
    Seff1.main=mat[,1]-covS[3,1]/covS[3,3]*mat[,3]
    #Seff1.interac=mat[,2]-covS[3,2]/covS[3,3]*mat[,3]
    Seff1.interac=mat[,2]-covS[1,2]/covS[1,1]*mat[,1]-covS[3,1]/covS[3,3]*mat[,3]
    
    #Seff1=mat[,1]+mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    
    Seff2=mat[,3]-covS[1,3]/covS[1,1]*mat[,1] #-covS[2,3]/covS[2,2]*mat[,2]
    
    if(sum(w)==0 | sum(w)==1){
      qzi_use=w[1]*Seff1.main+w[2]*Seff1.interac+(1-sum(w))*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a.main=round(n2samp*w[1])
      ka.main=round(n2a.main/2*0.15)
      res01.main <- rkd_ext_opt.f(ka.main, n2a.main, Seff1.main, phaseI_strat)
      
      n2a.interac=round(n2samp*w[2])
      ka.interac=round(n2a.interac/2*0.15)
      res01.interac <- rkd_ext_opt.f(ka.interac, n2a.interac, Seff1.interac, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% c(res01.main$s_id, res01.interac$s_id)), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a.main-n2a.interac
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01.main$s_m-res01.interac$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, c(res01.main$s_id, res01.interac$s_id))
      
    }
    
  }else if(design %in% "TRSD-v5-eff-interac-seq"){
    mat=as.matrix(cbind(qzi[[1]], qzi[[1]]*phaseI_strat$dt_ext$x1, qzi[[2]]))
    
    covS=cov(mat)
    Seff1.main=mat[,1]-covS[3,1]/covS[3,3]*mat[,3]
    #Seff1.interac=mat[,2]-covS[3,2]/covS[3,3]*mat[,3]
    Seff1.interac=mat[,2]-covS[1,2]/covS[1,1]*mat[,1] #-covS[3,1]/covS[3,3]*mat[,3]
    
    #Seff1=mat[,1]+mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    
    Seff2=mat[,3]-covS[1,3]/covS[1,1]*mat[,1] #-covS[2,3]/covS[2,2]*mat[,2]
    
    if(sum(w)==0 | sum(w)==1){
      qzi_use=w[1]*Seff1.main+w[2]*Seff1.interac+(1-sum(w))*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a.main=round(n2samp*w[1])
      ka.main=round(n2a.main/2*0.15)
      res01.main <- rkd_ext_opt.f(ka.main, n2a.main, Seff1.main, phaseI_strat)
      
      n2a.interac=round(n2samp*w[2])
      ka.interac=round(n2a.interac/2*0.15)
      res01.interac <- rkd_ext_opt.f(ka.interac, n2a.interac, Seff1.interac, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% c(res01.main$s_id, res01.interac$s_id)), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a.main-n2a.interac
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01.main$s_m-res01.interac$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, c(res01.main$s_id, res01.interac$s_id))
      
    }
    
  }else if(design %in% "FRSD-v5-eff-interac-seq"){
    mat=as.matrix(cbind(qzi[[1]], qzi[[1]]*phaseI_strat$dt_ext$x1, qzi[[2]], qzi[[2]]*phaseI_strat$dt_ext$x1))
    
    covS=cov(mat)
    Seff1.main=mat[,1]-covS[3,1]/covS[3,3]*mat[,3]
    #Seff1.interac=mat[,2]-covS[3,2]/covS[3,3]*mat[,3]
    Seff1.interac=mat[,2]-covS[1,2]/covS[1,1]*mat[,1] #-covS[3,1]/covS[3,3]*mat[,3]
    
    #Seff1=mat[,1]+mat[,2]-covS[3,1]/covS[3,3]*mat[,3]
    
    Seff2.main=mat[,3]-covS[1,3]/covS[1,1]*mat[,1] #-covS[2,3]/covS[2,2]*mat[,2]
    Seff2.interac=mat[,4]-covS[3,4]/covS[3,3]*mat[,3] #-covS[3,1]/covS[3,3]*mat[,3]
    
    if(sum(w[1:2])==0 | sum(w[1:2])==1){
      qzi_use=w[1]*Seff1.main+w[2]*Seff1.interac+w[3]*Seff2.main+w[4]*Seff2.interac
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a.main=round(n2samp*w[1])
      ka.main=round(n2a.main/2*0.15)
      res01.main <- rkd_ext_opt.f(ka.main, n2a.main, Seff1.main, phaseI_strat)
      
      n2a.interac=round(n2samp*w[2])
      ka.interac=round(n2a.interac/2*0.15)
      res01.interac <- rkd_ext_opt.f(ka.interac, n2a.interac, Seff1.interac, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% c(res01.main$s_id, res01.interac$s_id)), "R"] <- 0
      dt_b.main=subset(dt_a_use, R==0)
      
      n2b.main=round((n2samp-n2a.main-n2a.interac)*w[3])
      kb.main=round(n2b.main/2*0.15)
      phaseI_strat.b.main=phaseI_strat
      phaseI_strat.b.main$dt_ext=dt_b.main
      phaseI_strat.b.main$strata.n=phaseI_strat$strata.n-res01.main$s_m-res01.interac$s_m
      res02.main <- rkd_ext_opt_2st.f(kb.main, n2b.main, Seff2.main, phaseI_strat, c(res01.main$s_id, res01.interac$s_id))
      
      
      dt_a.b1_use <- phaseI_strat$dt_ext #dt_c
      dt_a.b1_use$R<-1
      dt_a.b1_use[!(dt_a_use$id %in% c(res01.main$s_id, res01.interac$s_id, res02.main$s_id)), "R"] <- 0
      dt_b.interac=subset(dt_a.b1_use, R==0)
      
      n2b.interac=n2samp-n2a.main-n2a.interac-n2b.main
      kb.interac=round(n2b.interac/2*0.15)
      phaseI_strat.b.interac=phaseI_strat
      phaseI_strat.b.interac$dt_ext=dt_b.interac
      phaseI_strat.b.interac$strata.n=phaseI_strat$strata.n-res01.main$s_m-res01.interac$s_m-res02.main$s_m
      
      res0 <- rkd_ext_opt_2st.f(kb.interac, n2b.interac, Seff2.interac, phaseI_strat, c(res01.main$s_id, res01.interac$s_id, res02.main$s_id))
    }
    
  }
  
  
  
  
  res = list(sel=res0, qzi=qzi)
  
  return(res)
}