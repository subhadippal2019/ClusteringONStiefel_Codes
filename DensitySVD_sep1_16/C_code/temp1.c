#include "stdio.h"
#include "stdlib.h"
#include "gsl/gsl_sf_exp.h"

//double gsl_sf_exp (double x)
void gslExp(double *x, double* result){
    /*
    This function returns the Fibonacci sequence
    where F[n] = F[n-1] + F[n-2]
    */
    double y ;
    y = gsl_sf_exp(*x);
    *result = y;

}
