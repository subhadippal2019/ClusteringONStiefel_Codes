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

//R CMD SHLIB del_hyper_2by2_R.c -I/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/include -L/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/lib -lgsl -lgslcblas -lm

//extern "C"
//{

	
	
	double calculate_A_k(double d2,double p,int k){

		double tmp1,tmp2,tmp3,tmp4;
		tmp1 = d2*d2/16;
		tmp2 = gsl_sf_poch(p-0.5,k)*gsl_sf_poch(p,2*k)*gsl_sf_fact(k);
		tmp3 = gsl_sf_gamma(p+2*k);

		tmp4 = tmp1*1.0/tmp2*tmp3;

		return tmp4;

	}

	

	double calculate_term1(double d1,double t_plus_c,double p,int k){

		double r,z,nu,T1,T2,T3,T4,T5;
		double a = p+(2.0*k);
		r = (1.0-a)*0.5;
		z = 2*sqrt(t_plus_c);
		nu = (a-1.0);

		T1 = r*(gsl_sf_log(t_plus_c));
		T2 = gsl_sf_bessel_Inu(nu,z);
		T3 = gsl_sf_log(T2);
		T4 = T1+T3;
		T5 = 2*d1*(gsl_sf_exp(T4));

		return T5;
	}
	/// term2 can be calculated from term1 easily

	double calculate_term3b(double d1,double t_plus_c,double p,int k){

		double r,z,nu,T1,T2,T3,T4,T5;
		double a = p+(2.0*k);
		r = -0.5*a;
		z = 2*sqrt(t_plus_c);
		nu = a;


		T1 = r*(gsl_sf_log(t_plus_c));
		T2 = gsl_sf_bessel_Inu(nu,z);
		T3 = gsl_sf_log(T2);
		T4 = T1+T3;
		T5 = d1*d1*d1*0.5*(gsl_sf_exp(T4));

		return T5;
	}


	

	double calculate_del_hyper_2by2(double p,double d1,double d2,int KK){

		double S=0.0;
		int k=0;

		double c,t,t_plus_c,A_k,a;
		double term1,term2,term3a,term3b,sum_T;

		c = d2*d2/4;
		t = d1*d1/4;
		t_plus_c = t+c;
		

		for(k=0;k<KK;k++){

			term1 = calculate_term1(d1,t_plus_c,p,k);
		  //a = p+(2.0*k);
			//term2 = d1*d1*0.25*term1*(1.0/(t+c));
			//term3a = term2*(a-1.0)*0.5;
			term3b = calculate_term3b(d1,t_plus_c,p,k);
			A_k = calculate_A_k(d2,p,k);
			sum_T = term1+term3b;

			S += A_k*sum_T;

		}

		return S;
	}

	void del_hyper_2by2_R(double *cc,double *D, double *dRet){

		int KK = 30; // what to choose ??
		
		double p = *cc;
		double d1 = 2*sqrt(D[0]);
		double d2 = 2*sqrt(D[1]);
		

		*dRet = calculate_del_hyper_2by2(p,d1,d2,KK);	

		
	}
		
//}
