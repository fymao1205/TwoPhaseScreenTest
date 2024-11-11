
# max. comp loglik



# max. obs loglik
est.obs_pwc.f <- function(dat, interac.ind=c(0,0), brks){

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
  
  obj.f <- function(para){
    
    # alp = exp(para[1:len]) #c(exp(para[1]), 1)
    # beta = para[len+(1:2)]
    # eta = para[len+2+(1:3)]
    # gam = para[len+6]
    
    #eta0=exp(para[1:3])
    #eta = eta0/(1+eta0) #para[1:2]
    eta = para[1:2]
    p1 = expit.f(para[3])
    alp = exp(para[3+1:len])
    beta = c(para[3+len+1:len.beta], 0)
    gam = c(para[3+len+len.beta+1:len.gam],0)
    
    #res0 <- logLn_obs_pwc(dat$r, dat$rz, dat$z.o, dat$d, dat$a, dat$x1, dat$x2, gam, alp, beta, eta, brks)
    #res0 <- logLn_obs_pwc(dat$r, dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, eta, p1, brks)
    res0 <- RtoC_logLn_obs_pwc(dat$r, dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, eta, p1, brks)

    res <- -sum(res0)
    
    # print(para)
    # print(res)
    
    return(res)
    
  }
  
  est <- optim(par=rep(0.01,len+3+len.beta+len.gam), obj.f, method=("L-BFGS-B"), hessian=TRUE)
  #est <- nlm(f=obj.f, p=rep(0.01, p1+p2+p3+4), hessian=TRUE)
  return(est)
}


logL_pwc.f <- function(para, dat, interac.ind=c(0,0), brks){
  
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

  #res0 <- logLn_obs_pwc(dat$r, dat$rz, dat$z.o, dat$d, dat$a, dat$x1, dat$x2, gam, alp, beta, eta, brks)
  #res0 <- logLn_obs_pwc(dat$r, dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, eta, p1, brks)
  res0 <- RtoC_logLn_obs_pwc(dat$r, dat$d1, dat$d2, dat$a1, dat$a2, dat$x1, dat$x2, gam, alp, beta, eta, p1, brks)
  
  return(res0)
}

score_pwc.f <- function(grad=1e-06, para, dat, interac.ind=c(0,0), brks){
  
  logL0 <- logL_pwc.f(para, dat, interac.ind, brks)
  sc.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    logL.p <- logL_pwc.f(para.p, dat, interac.ind, brks)
    (logL0-logL.p)/grad 
  })
  
  res <- do.call("cbind", sc.list)
  return(res)
}

hess_pwc.f <- function(grad=1e-06, para, dat, interac.ind, brks){
  
  s0 <- score_pwc.f(grad, para, dat, interac.ind, brks)
  h.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    s.p <- score_pwc.f(grad, para.p, dat, interac.ind, brks)
    (colSums(s.p)-colSums(s0))/grad 
  })
  
  res <- do.call("cbind", h.list)
  return(res)
  
}

ase_pwc.f <- function(grad=1e-06, est, dat, interac.ind, brks){
  # score vector
  score.mat = score_pwc.f(grad, est$par, dat, interac.ind, brks)
  #score.mat = score_trueh.f(grad, est$par, del, time, v.mat, pos.cov.T, pos.cov.Z, pos.cov.X)
  
  #hess.mat = hess_pwc.f(grad, est, dat, interac.ind, brks)
  hess.mat = est$hess 
  #hessian(lik.f, est, x1=dt.cov[,1], v=dt.cov[,-1, drop=FALSE], del1, del2, b1, b2, cov_dist) # hessian matrix by numeircal approximation
  B = t(score.mat) %*% score.mat
  A = (hess.mat)
  #A = (hesmat)
  avar = ginv(A)%*% B%*% ginv(A)
  ase = sqrt(diag(avar))
  res = list(ase=ase, avar=avar, A=A, B=B, score=colSums(score.mat))
  return(res)
}

