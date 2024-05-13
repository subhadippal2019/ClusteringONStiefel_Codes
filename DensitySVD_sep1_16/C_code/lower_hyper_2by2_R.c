#include "stdio.h"
#include "R.h"
#include "Rmath.h"
#include "stdlib.h"
#include "gsl/gsl_sf_pow_int.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_hyperg.h"

//double gsl_sf_pow_int (double x, int n)
//double gsl_sf_fact (unsigned int n)
//double gsl_sf_poch (double a, double x)
//double gsl_sf_hyperg_0F1 (double c, double x)

// R CMD SHLIB lower_hyper_2by2_R.c -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm
// R CMD SHLIB lower_hyper_2by2_R.c -I/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/include -L/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/lib -lgsl -lgslcblas -lm

//extern "C"
//{

	void lower_hyper_2by2_R(double *cc,double *R,double *d){

		int KK = 30;
		double S=0.0;
		double tmp1,tmp2,tmp3,tmp4;
		int k=0;
		
		double r1=R[0];
		double r2=R[1];
		
		double c = *cc;

		for(k=0;k<KK;k++){
		
			tmp1 = r1*r2*1.0;

			tmp2 = gsl_sf_fact(k);
			tmp3 = gsl_sf_fact(2*k+1);

			S += gsl_sf_pow_int(tmp1,k)/(tmp2*tmp2*tmp3);

		}

		*d = S;
	}

	void lower_0F1_scalar(double *a,double *z, double*d){
	
		double x;
		x = gsl_sf_hyperg_0F1(*a,*z);
		*d = x;
	}


