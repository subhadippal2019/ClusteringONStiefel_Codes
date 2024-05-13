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

//R CMD SHLIB hyper_2by2_R.c -I/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/include -L/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/lib -lgsl -lgslcblas -lm

//extern "C"
//{

	void partial_hyper_2by2_R(double *c,double *D,double *d){

		int KK = 30;
		double S=0.0;
		double tmp1,tmp2,tmp3,tmp4,tmp5;
		int k=0;
		
		double r1=D[0];
		double r2=D[1];
		
		double d1 = 2*sqrt(r1);
		double d2 = 2*sqrt(r2);


		double p = *c;

		for(k=0;k<KK;k++){
		
			tmp1 = r1*r2*1.0;

			tmp2 = gsl_sf_poch(p-0.5,k);
			tmp3 = gsl_sf_poch(p,2*k);
			tmp4 = gsl_sf_fact(k);

			tmp5 = tmp2*tmp3*tmp4;		

			S += gsl_sf_pow_int(tmp1,k)*(   ((2*k/d1)*gsl_sf_hyperg_0F1(p+2*k,r1+r2)) + 
										    ((d1/(2*p+4*k))*gsl_sf_hyperg_0F1(p+2*k+1,r1+r2))   )*1.0/tmp5;

		}

		*d = S;
	}

	void partial_partial_hyper_2by2_R(double *c,double *D,double *d){

		int KK = 30;
		double S=0.0;
		double tmp1,tmp2,tmp3,tmp4,tmp5;
		int k=0;
		
		double r1=D[0];
		double r2=D[1];
		
		double d1 = 2*sqrt(r1);
		double d2 = 2*sqrt(r2);


		double p = *c;

		for(k=0;k<KK;k++){
		
			tmp1 = r1*r2*1.0;

			tmp2 = gsl_sf_poch(p-0.5,k);
			tmp3 = gsl_sf_poch(p,2*k);
			tmp4 = gsl_sf_fact(k);

			tmp5 = tmp2*tmp3*tmp4;		

			S += gsl_sf_pow_int(tmp1,k)*(  ((2*k)*(2*k-1)/(d1*d1)*gsl_sf_hyperg_0F1(p+2*k,r1+r2)) + 
										   ((4*k+1)/(2*p+4*k) *gsl_sf_hyperg_0F1(p+2*k+1,r1+r2)) +  
										   ((d1*d1)/((2*p+4*k)*(2*p+4*k+2))*gsl_sf_hyperg_0F1(p+2*k+2,r1+r2))  )*1.0/tmp5;

		}

		*d = S;
	}


//}
