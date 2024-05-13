//#include "stdio.h"
//#include "stdlib.h"
#include "mex.h"
#include "gsl/gsl_sf_pow_int.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_hyperg.h"

//double gsl_sf_pow_int (double x, int n)
//double gsl_sf_fact (unsigned int n)
//double gsl_sf_poch (double a, double x)
//double gsl_sf_hyperg_0F1 (double c, double x)


//mex hyper_2by2_MEX.c -I/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/include  -L/Users/subhajit/SOFTWARES/gsl-2.2.1/GSL_install/lib -lgsl -lgslcblas -lm

//double calculate_hyper_2by2(double c,double r1,double r2,int KK);
/*
int main(int argc,char* argv[]){

	int K=0;
	double r1=0.0,r2=0.0,d=0.0,c=0.0;
	if(argc < 3){
		fprintf(stderr,"\nUsage: ./hyper_2_by_2 <diag 1> <diag 2>\n\n");
		exit(0);
	}

	
	sscanf(argv[1],"%lf",&r1);
	sscanf(argv[2],"%lf",&r2);
	
	K = 30;
	c = 1.5;
	d = calculate_hyper_2by2(c,r1,r2,K);

	fprintf(stdout,"0F1(3/2;[%lf %lf]) = %lf\n\n",r1,r2,d);

	return 0;
}
*/

void hyper_2by2_MEX(double c,double* R,double* d,mwSize KK){

	mwSize k;
	double S=0.0;
	double tmp1,tmp2,tmp3,tmp4,tmp5;

	double r1 = R[0];
	double r2 = R[1];
	
	for(k=0;k<KK;k++){
		
		tmp1 = r1*r2*1.0;

		tmp2 = gsl_sf_poch(c-0.5,k);
		tmp3 = gsl_sf_poch(c,2*k);
		tmp4 = gsl_sf_fact(k);

		tmp5 = tmp2*tmp3*tmp4;		

		S += gsl_sf_pow_int(tmp1,k)*gsl_sf_hyperg_0F1(c+2*k,r1+r2)*1.0/tmp5;

	}

	*d = S;
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

	if(nrhs != 2) {
    	mexErrMsgIdAndTxt("Hyper_2by2:nrhs",
							 "Two inputs required.");
	}

	if(nlhs != 1) {
    	mexErrMsgIdAndTxt("Hyper_2by2:nlhs",
                      		"One output required.");
	}

	if( !mxIsDouble(prhs[1])){ 
    	mexErrMsgIdAndTxt("Hyper_2by2:notDouble",
        			"Input matrix must be type double.");
	}
	
	if(mxGetM(prhs[1]) != 1) {
    	mexErrMsgIdAndTxt("Hyper_2by2:notRowVector",
            			    "Input must be a row vector.");
	}

	double c;      /* input scalar */
	double *R;       /* 1xN input matrix */

	mwSize K;

	double *d;

	c = mxGetScalar(prhs[0]);
	R = mxGetPr(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

	d = mxGetPr(plhs[0]);

	K = 30;

	hyper_2by2_MEX(c,R,d,K);

}


