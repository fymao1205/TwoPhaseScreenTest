//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>

#include "commonf.h"

using namespace Rcpp;
using namespace arma;


// Obtain environment containing function
//Rcpp::Environment package_env("package:mnormt"); 

// Make function callable from C++
//Rcpp::Function bipmvnrm = package_env["pmnorm"];    

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export()]]
double logL_csd_cal_wei(int rz, int z, int d, double t, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  double Hf = pow(alp[0]*t, alp[1])*exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2);
  double Sf = exp(-Hf); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(rz==1){
    
    if(d==1){
      logL1=z*log(piZ) + (1-z)*(log(1-Sf)+log(1-piZ));
    }else{
      logL1=(-Hf+log(1-piZ));
    }
    //logL1 = d*z*log(piZ) + d*(1-z)*(log(1-Sf)+log(1-piZ)) + (1-d)*(-Hf+log(1-piZ));
  }else{
    //d*(p.z1+(1-surv.a)*p.z0)+(1-d)*surv.a*p.z0
    logL1 = d*log(piZ+(1-Sf)*(1-piZ))+(1-d)*(-Hf+log(1-piZ));
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  double logL2 = x1*(x2*(eta[0]+eta[1]) - log(1+exp(eta[0]+eta[1])) + log(p1) ) + (1-x1)*( x2*(eta[0]) - log(1+exp(eta[0])) + log(1-p1) );
  
  double res=logL1+logL2;
  
  return res;
}

// [[Rcpp::export()]]
double logL_csd_obs_wei(int r, int rz, int z, int d, double t, double x1, double x2, 
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double res=0.0;
  if(r==1){
    res = logL_csd_cal_wei(rz, z, d, t, x1, x2, gam, alp, beta, eta, p1);
  }else{
    double res0 = logL_csd_cal_wei(rz, z, d, t, x1, 0.0, gam, alp, beta, eta, p1);
    double res1 = logL_csd_cal_wei(rz, z, d, t, x1, 1.0, gam, alp, beta, eta, p1);
    res = log(exp(res0)+exp(res1));
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_csd_cal_wei(IntegerVector rz, IntegerVector z, IntegerVector d, NumericVector t, NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, 
                            NumericVector eta, double p1){
  
  NumericVector res=no_init(d.size());
  
  for(int i=0; i<d.size(); i++){
    res[i] = logL_csd_cal_wei(rz[i], z[i], d[i], t[i], x1[i], x2[i], gam, alp, beta, eta, p1); //wgt[i]*(temp);
  }
  
  return res;
}

// [[Rcpp::export()]]
double logL_cal_wei(int d1, int d2, double t1, double t2, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  double expcovT=exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2);
  double Hft1 = pow(alp[0]*t1, alp[1])*expcovT;
  double Sft1 = exp(-Hft1); 
  double Hft2 = pow(alp[0]*t2, alp[1])*expcovT;
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        //logL1=(log(1-piZ)-Hft2);
        logL1=0;
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  double logL2 = x1*(x2*(eta[0]+eta[1]) - log(1+exp(eta[0]+eta[1])) + log(p1) ) + (1-x1)*( x2*(eta[0]) - log(1+exp(eta[0])) + log(1-p1) ); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  //double logL2 = x1*x2*log(eta[0]) + x1*(1-x2)*log(eta[1]) + (1-x1)*x2*log(eta[2]) + (1-x1)*(1-x2)*log(1-eta[0]-eta[1]-eta[2]);
  
  double res=logL1+logL2;
  
  return res;
}

// [[Rcpp::export()]]
double logL_cal1_wei(int d1, int d2, double t1, double t2, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  double expcovT=exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2);
  double Hft1 = pow(alp[0]*t1, alp[1])*expcovT;
  double Sft1 = exp(-Hft1); 
  double Hft2 = pow(alp[0]*t2, alp[1])*expcovT;
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        logL1=0;//(log(1-piZ)-Hft2);
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  
  double res=logL1;
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_cal1_wei(IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, 
                             NumericVector x1, NumericVector x2, 
                             NumericVector gam, NumericVector alp, NumericVector beta){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_cal1_wei(d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta); //wgt[i]*(temp);
  }
  
  return res;
}

extern "C" {

  void C_logLn_cal1_wei(int *R_d1, int *R_d2, double *R_t1, double *R_t2, 
                        double *R_x1, double *R_x2, 
                        double *R_gam, double *R_alp, double *R_beta,
                        int *R_d1_size, int *R_gam_size, int *R_alp_size, int *R_beta_size,
                        double *ret) {

    int n      = *R_d1_size;

    IntegerVector d1(n);
    IntegerVector d2(n);
    NumericVector t1(n);
    NumericVector t2(n);
    NumericVector x1(n);
    NumericVector x2(n);
    NumericVector gam(*R_gam_size);
    NumericVector alp(*R_alp_size);
    NumericVector beta(*R_beta_size);

    copyIntToIntVec(R_d1, n, d1);
    copyIntToIntVec(R_d2, n, d2);
    copyDblToNumVec(R_t1, n, t1);
    copyDblToNumVec(R_t2, n, t2);
    copyDblToNumVec(R_x1, n, x1);
    copyDblToNumVec(R_x2, n, x2);
    copyDblToNumVec(R_gam, *R_gam_size, gam);
    copyDblToNumVec(R_alp, *R_alp_size, alp);
    copyDblToNumVec(R_beta, *R_beta_size, beta);

    NumericVector vec = logLn_cal1_wei(d1, d2, t1, t2, x1, x2, gam, alp, beta);
    for (int i=0; i<n; i++) ret[i] = vec[i];

    return;
  }
}

double logL_obs_wei(int r, int d1, int d2, double t1, double t2, double x1, double x2, 
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  double res=0.0;
  if(r==1){
    res = logL_cal_wei(d1, d2, t1, t2, x1, x2, gam, alp, beta, eta, p1);
  }else{
    double res0 = logL_cal_wei(d1, d2, t1, t2, x1, 0.0, gam, alp, beta, eta, p1);
    double res1 = logL_cal_wei(d1, d2, t1, t2, x1, 1.0, gam, alp, beta, eta, p1);
    res = log(exp(res0)+exp(res1));
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_obs_wei(IntegerVector r, IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_obs_wei(r[i], d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, eta, p1); //wgt[i]*(temp);
  }
  
  return res;
}

// [[Rcpp::export()]]
double logL_cal_pwc(double d1, double d2, double t1, double t2, double x1, double x2,
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1,  
                    NumericVector brks){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); 
  NumericVector levs = alp*exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2); 
  double Hft1 = Hpwc_double(t1, brks, levs, 0);
  double Sft1 = exp(-Hft1); 
  double Hft2 = Hpwc_double(t2, brks, levs, 0);
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        logL1=0;//(log(1-piZ)-Hft2);
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  //double logL2 = x1*x2*log(eta[0]) + x1*(1-x2)*log(eta[1]) + (1-x1)*x2*log(eta[2]) + (1-x1)*(1-x2)*log(1-eta[0]-eta[1]-eta[2]);
  double logL2 = x1*(x2*(eta[0]+eta[1]) - log(1+exp(eta[0]+eta[1])) + log(p1) ) + (1-x1)*( x2*eta[0] - log(1+exp(eta[0])) + log(1-p1) );
  
  double res=logL1+logL2;
  
  return res;
}

extern "C" {

  void C_logL_cal_pwc(double *R_d1, double *R_d2, double *R_t1, double *R_t2, double *R_x1, double *R_x2,
                      double *R_gam, double *R_alp, double *R_beta, double *R_eta, double *R_p1,  
                      double *R_brks, int *R_gam_size, int *R_alp_size, int *R_beta_size,
                      int *R_eta_size, int *R_brks_size, double *ret) {

    NumericVector gam(*R_gam_size);
    NumericVector alp(*R_alp_size);
    NumericVector beta(*R_beta_size);
    NumericVector eta(*R_eta_size);
    NumericVector brks(*R_brks_size);

    copyDblToNumVec(R_gam,  *R_gam_size,  gam);
    copyDblToNumVec(R_alp,  *R_alp_size,  alp);
    copyDblToNumVec(R_beta, *R_beta_size, beta);
    copyDblToNumVec(R_eta,  *R_eta_size,  eta);
    copyDblToNumVec(R_brks, *R_brks_size, brks);

    *ret = logL_cal_pwc(*R_d1, *R_d2, *R_t1, *R_t2, *R_x1, *R_x2,
                        gam, alp, beta, eta, *R_p1,  brks);
    return;
  }

}

// [[Rcpp::export()]]
double logL_cal1_pwc(int d1, int d2, double t1, double t2, double x1, double x2,
                     NumericVector gam, NumericVector alp, NumericVector beta, 
                     NumericVector brks){
  double expcovZ=exp(gam[0]+x1*gam[1]+x2*gam[2]+gam[3]*x1*x2); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  //double Hf = pow(alp[0]*t, alp[1])*exp(beta[0]*x1+beta[1]*x2);
  NumericVector levs = alp*exp(beta[0]*x1+beta[1]*x2+beta[2]*x1*x2); 
  //double hf = hpwc_double(t, brks, levs, 0); 
  double Hft1 = Hpwc_double(t1, brks, levs, 0);
  double Sft1 = exp(-Hft1); 
  double Hft2 = Hpwc_double(t2, brks, levs, 0);
  double Sft2 = exp(-Hft2); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1;
  if(d1==1){
    
    logL1=log(piZ);
  }else{
    
    if(d1>1){
      //logL1=d2*log((1-Sft2)*(1-piZ)+piZ)+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log((1-Sft2)*(1-piZ)+piZ);
      }else{
        logL1=0;//(log(1-piZ)-Hft2);
      }
    }else{//d1==0 or 1
      //logL1=d2*(log( Sft1-Sft2 )+log(1-piZ))+(1-d2)*(log(1-piZ)-Hft2);
      if(d2==1){
        logL1=log( Sft1-Sft2 )+log(1-piZ);
      }else{
        logL1=(log(1-piZ)-Hft2);
      }
    }
  }
  
  //double logL2 = x2*(eta[0]+eta[1]*x1) - log(1+exp(eta[0]+eta[1]*x1)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  
  double res=logL1;//+logL2;
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_cal1_pwc(IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, 
                            NumericVector brks){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_cal1_pwc(d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, brks); //wgt[i]*(temp);
  }
  
  return res;
}


// [[Rcpp::export()]]
NumericVector logLn_cal_pwc(IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, NumericVector x1, NumericVector x2, 
                             NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1,  
                             NumericVector brks){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_cal_pwc(d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, eta, p1, brks); //wgt[i]*(temp);
  }
  
  return res;
}


double logL_obs_pwc(int r, int d1, int d2, double t1, double t2, double x1, double x2, 
                    NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1, 
                    NumericVector brks){
  double res=0.0;
  if(r==1){
    res = logL_cal_pwc(d1, d2, t1, t2, x1, x2, gam, alp, beta, eta, p1, brks);
  }else{
    double res0 = logL_cal_pwc(d1, d2, t1, t2, x1, 0.0, gam, alp, beta, eta, p1, brks);
    double res1 = logL_cal_pwc(d1, d2, t1, t2, x1, 1.0, gam, alp, beta, eta, p1, brks);
    res = log(exp(res0)+exp(res1));
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_obs_pwc(IntegerVector r, IntegerVector d1, IntegerVector d2, NumericVector t1, NumericVector t2, 
                            NumericVector x1, NumericVector x2, 
                            NumericVector gam, NumericVector alp, NumericVector beta, NumericVector eta, double p1, 
                            NumericVector brks){
  
  NumericVector res=no_init(d1.size());
  
  for(int i=0; i<d1.size(); i++){
    res[i] = logL_obs_pwc(r[i], d1[i], d2[i], t1[i], t2[i], x1[i], x2[i], gam, alp, beta, eta, p1, brks); //wgt[i]*(temp);
  }
  
  return res;
}

extern "C" {

  void C_logLn_obs_pwc(int *R_r, int *R_d1, int *R_d2, double *R_t1, double *R_t2, 
                       double *R_x1, double *R_x2, double *R_gam, double *R_alp, 
                       double *R_beta, double *R_eta, double *R_p1, double *R_brks,
                       int *R_d1_size, int *R_gam_size, int *R_alp_size, int *R_beta_size,
                       int *R_eta_size, int *R_brks_size, double *ret) {

    int n = *R_d1_size;
    double p1 = *R_p1;

    IntegerVector r(n);
    IntegerVector d1(n);
    IntegerVector d2(n);
    NumericVector t1(n);
    NumericVector t2(n);
    NumericVector x1(n);
    NumericVector x2(n);
    NumericVector gam(*R_gam_size);
    NumericVector alp(*R_alp_size);
    NumericVector beta(*R_beta_size);
    NumericVector eta(*R_eta_size);
    NumericVector brks(*R_brks_size);

    copyIntToIntVec(R_r, n, r);
    copyIntToIntVec(R_d1, n, d1);
    copyIntToIntVec(R_d2, n, d2);
    copyDblToNumVec(R_t1, n, t1);
    copyDblToNumVec(R_t2, n, t2);
    copyDblToNumVec(R_x1, n, x1);
    copyDblToNumVec(R_x2, n, x2);
    copyDblToNumVec(R_gam, *R_gam_size, gam);
    copyDblToNumVec(R_alp, *R_alp_size, alp);
    copyDblToNumVec(R_beta, *R_beta_size, beta);
    copyDblToNumVec(R_eta, *R_eta_size, eta);
    copyDblToNumVec(R_brks, *R_brks_size, brks);

    NumericVector vec = logLn_obs_pwc(r, d1, d2, t1, t2, x1, x2, gam, alp, beta, eta, p1, brks);
    for (int i=0; i<n; i++) ret[i] = vec[i];

    return;
  }

}

// [[Rcpp::export()]]
NumericVector vHpwc_H1(NumericVector t, NumericVector x1, NumericVector x2, 
                       NumericVector alp, NumericVector beta, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(beta[0]*x1[i]+beta[1]*x2[i]);
    res[i] = Hpwc_double(t[i], brks, levs, 0); 
  }
  
  return res;
  
}

// [[Rcpp::export()]]
NumericVector vhpwc_H1(NumericVector t, NumericVector x1, NumericVector x2, 
                       NumericVector alp, NumericVector beta, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(beta[0]*x1[i]+beta[1]*x2[i]);
    res[i] = hpwc_double(t[i], brks, levs, 0); 
  }
  
  return res;
  
}


// [[Rcpp::export()]]
NumericVector vppwc_H1(NumericVector t, NumericVector x1, NumericVector x2, 
                       NumericVector alp, NumericVector beta, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(beta[0]*x1[i]+beta[1]*x2[i]);
    res[i] = ppwc_double(t[i], brks, levs, 1,0); 
  }
  
  return res;
  
}

