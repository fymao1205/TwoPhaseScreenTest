twoPhaseScreen_predict <- function(est.obj, t=0, grad=1e-06) {

  check_est.obj(est.obj) 
  check_num(grad, "grad")

  ret0 <- ase.cumrisk_g_x2_pwc.f(grad=grad, para.varmat=est.obj$cov, t=t, x2=0, 
            para=est.obj$est, est.obj$brks, interac.ind=est.obj$interac.ind)
  ret1 <- ase.cumrisk_g_x2_pwc.f(grad=grad, para.varmat=est.obj$cov, t=t, x2=1, 
            para=est.obj$est, est.obj$brks, interac.ind=est.obj$interac.ind)

  list(est.cNPV=ret0$est, se.cNPV=ret0$ase, est.PPV=ret1$est, se.PPV=ret1$ase)
}
