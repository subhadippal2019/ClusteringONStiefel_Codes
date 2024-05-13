#include "stdio.h"
#include "stdlib.h"

#include "gsl/gsl_sf_pow_int.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_hyperg.h"

//double gsl_sf_pow_int (double x, int n)
//double gsl_sf_fact (unsigned int n)
//double gsl_sf_poch (double a, double x)
//double gsl_sf_hyperg_0F1 (double c, double x)

//gcc -L/usr/local/lib example.o -lgsl -lgslcblas -lm

double calculate_hyper_2by2(double c,double r1,double r2,int KK);

int main(int argc,char* argv[]){

	int K=0;
	double r1=0.0,r2=0.0,d=0.0,c=0.0;
	if(argc < 4){
		fprintf(stderr,"\nUsage: ./hyper_2_by_2 <c> <diag 1> <diag 2>\n\n");
		exit(0);
	}

	
	sscanf(argv[1],"%lf",&c);
	sscanf(argv[2],"%lf",&r1);
	sscanf(argv[3],"%lf",&r2);
	
	K = 30;
	
	d = calculate_hyper_2by2(c,r1,r2,K);

	fprintf(stdout,"0F1(3/2;[%lf %lf]) = %lf\n\n",r1,r2,d);

	return 0;
}

double calculate_hyper_2by2(double c,double r1,double r2,int KK){

	double S=0.0;
	double tmp1,tmp2,tmp3,tmp4,tmp5;
	int k=0;
	tmp1 = r1*r2*1.0;
	for(k=0;k<KK;k++){

		tmp2 = gsl_sf_poch(c-0.5,k);
		tmp3 = gsl_sf_poch(c,2*k);
		tmp4 = gsl_sf_fact(k);

		tmp5 = tmp2*tmp3*tmp4;		

		S += gsl_sf_pow_int(tmp1,k)*gsl_sf_hyperg_0F1(c+2*k,r1+r2)*1.0/tmp5;

	}

	return S;
}

