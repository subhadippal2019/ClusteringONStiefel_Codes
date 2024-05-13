#include "stdio.h"
#include "R.h"
#include "Rmath.h"
#include "stdlib.h"
#include "gsl/gsl_sf_pow_int.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_log.h"
#include "gsl/gsl_sf_exp.h"
//double gsl_sf_pow_int (double x, int n)
//double gsl_sf_fact (unsigned int n)
//double gsl_sf_poch (double a, double x)
//double gsl_sf_hyperg_0F1 (double c, double x)
//double gsl_sf_bessel_Inu (double nu, double x)
//double gsl_sf_gamma (double x)
//double gsl_sf_log (double x)
//double gsl_sf_exp (double x)

//gcc -L/usr/local/lib example.o -lgsl -lgslcblas -lm

// this function calculated del(0F1(a;D^2/4))/del d1

//R CMD SHLIB del_hyper_2by2_R_alt.c -I/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/include -L/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/lib -lgsl -lgslcblas -lm

//extern "C"
//{
  
  
  
  double calculate_A_k(double d2,double p,int k){
    
    double T1,T2,T3,T4,A1_k;
    T1 = gsl_sf_pow_int(0.5*d2,2*k);
  	A1_k = 1.0/(gsl_sf_poch(p-0.5,k)*gsl_sf_poch(p,2*k)*gsl_sf_fact(k));
    T2 = gsl_sf_gamma(p+2*k);
    T3 = gsl_sf_pow_int(2.0,(p-1));
    
    T4 = T1*A1_k*T2*T3;
    
    return T4;
    
  }
  
  
  
  double calculate_term1(double d1,double d2,double p,int k){
    
    double g_d1,nu,T1,T2,T3_0,T3,T4,T5;

    g_d1 = sqrt(d1*d1+d2*d2);
    nu = (p+2*k-1.0);
    
    T1 = (2*k-1)*(gsl_sf_log(d1));
    T2 = -1.0*nu*(gsl_sf_log(g_d1));
    T3_0 = gsl_sf_bessel_Inu(nu,g_d1);
    T3 = gsl_sf_log(T3_0);

    T4 = T1+T2+T3;
    T5 = 2*k*(gsl_sf_exp(T4));
    
    return T5;
  }
  
  double calculate_term2(double d1,double d2,double p,int k){
    
    double g_d1,nu,nu_plus_1,T1,T2,T3_0,T3,T4,T5,T6;

    g_d1 = sqrt(d1*d1+d2*d2);
    nu = (p+2*k-1.0);
    
    T1 = (2*k)*(gsl_sf_log(d1));
    T2 = -1.0*nu*(gsl_sf_log(g_d1));
    nu_plus_1 = (nu+1.0);
    T3_0 = gsl_sf_bessel_Inu(nu_plus_1,g_d1);
    T3 = gsl_sf_log(T3_0);
    T4 = gsl_sf_log(d1) - gsl_sf_log(g_d1);

    T5 = T1+T2+T3+T4;
    T6 = gsl_sf_exp(T5);
    
    return T6;
  }
  
  
  
  
  double calculate_del_hyper_2by2(double p,double d1,double d2,int KK){
    
    double S=0.0;
    int k=0;
    
    double T1,T2,A_k;
    
    for(k=0;k<KK;k++){

    	A_k = calculate_A_k(d2,p,k);
    	T1 = calculate_term1(d1,d2,p,k);
    	T2 = calculate_term2(d1,d2,p,k);

    	S += A_k*(T1+T2);
    }
    
    return S;
  }
  
  void del_hyper_2by2_R_alt(double *cc,double *D, double *dRet){
    
    int KK = 30; // what to choose ??
    
    double p = *cc;
    double d1 = 2*sqrt(D[0]);
    double d2 = 2*sqrt(D[1]);
    
    
    *dRet = calculate_del_hyper_2by2(p,d1,d2,KK);	
    
    
  }
  
  //}
