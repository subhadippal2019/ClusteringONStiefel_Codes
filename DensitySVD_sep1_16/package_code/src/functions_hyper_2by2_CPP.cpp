

#include <stdio.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double hyper_2by2_R(double c,NumericVector R, double MaxTolarableError=0.0000000000001){
  int KK = 30;
  double S=0.0;
  //NumericVector term(KK);
  double temp,tmp1,tmp2,tmp3,tmp4,tmp5,num1,num2;
  int k=0;
  
  double r1=R[0];
  double r2=R[1];
  Rcpp::Environment gsl_env("package:gsl");
  
  Rcpp::Function gsl_sf_hyperg_0F1 = gsl_env["hyperg_0F1"];
  Rcpp::Function gsl_sf_fact= gsl_env["fact"];
  Rcpp::Function gsl_sf_poch= gsl_env["poch"];
  Rcpp::Function gsl_sf_pow_int= gsl_env["pow_int"];
  
  
  for(k=0;k<KK;k++){
    
    tmp1 = r1*r2*1.0;
    
    tmp2 = as<double> (gsl_sf_poch(c-0.5,k));
    tmp3 = as<double> ( gsl_sf_poch(c,2*k));
    tmp4 = as<double> (gsl_sf_fact(k));
    
    tmp5 = tmp2*tmp3*tmp4;		
    num1=  as<double> (gsl_sf_pow_int(tmp1,k));
    num2=  as<double> (gsl_sf_hyperg_0F1(c+2*k,r1+r2));
    //term[k]=(num1*num2*1.0/tmp5);
    temp=(num1*num2*1.0/tmp5);
    S += temp;
    if(  temp< MaxTolarableError) {
      /* terminate the loop using break statement */
        break;
    }
    
  }
  
  return S;
}


// [[Rcpp::export]]
NumericVector hyper_2by2_R_test(double c,NumericVector R, double MaxTolarableError=0.0000000000001){
  //double MaxTolarableError==0.0000000000001;
  int KK = 30;
  double S=0.0;
  NumericVector term(KK);
  double tmp1,tmp2,tmp3,tmp4,tmp5,num1,num2;
  int k=0;
  
  double r1=R[0];
  double r2=R[1];
  Rcpp::Environment gsl_env("package:gsl");
  
  Rcpp::Function gsl_sf_hyperg_0F1 = gsl_env["hyperg_0F1"];
  Rcpp::Function gsl_sf_fact= gsl_env["fact"];
  Rcpp::Function gsl_sf_poch= gsl_env["poch"];
  Rcpp::Function gsl_sf_pow_int= gsl_env["pow_int"];
  
  
  for(k=0;k<KK;k++){
    
    tmp1 = r1*r2*1.0;
    
    tmp2 = as<double> (gsl_sf_poch(c-0.5,k));
    tmp3 = as<double> ( gsl_sf_poch(c,2*k));
    tmp4 = as<double> (gsl_sf_fact(k));
    
    tmp5 = tmp2*tmp3*tmp4;		
    num1= as<double> (gsl_sf_pow_int(tmp1,k));
    num2=  as<double> (gsl_sf_hyperg_0F1(c+2*k,r1+r2));
    term[k]=(num1*num2*1.0/tmp5);
    S += (num1*num2*1.0/tmp5);
    if(  term[k]< MaxTolarableError) {
      /* terminate the loop using break statement */
        break;
    }
  }
  
  return term;
}



// [[Rcpp::export]]
double partial_d1_hyper_2by2_R(double c, NumericVector D, double MaxTolarableError=0.0000000000001){
  
  int KK = 30;
  double S=0.0;
  double temp,tmp1,tmp2,tmp3,tmp4,tmp5,num1,num2, num3;
  int k=0;
  
  double r1=D[0];
  double r2=D[1];
  
  Rcpp::Environment gsl_env("package:gsl");
  Rcpp::Function gsl_sf_hyperg_0F1 = gsl_env["hyperg_0F1"];
  Rcpp::Function gsl_sf_fact= gsl_env["fact"];
  Rcpp::Function gsl_sf_poch= gsl_env["poch"];
  Rcpp::Function gsl_sf_pow_int= gsl_env["pow_int"];
  
  double d1 = 2*sqrt(r1);
  //double d2 = 2*sqrt(r2);
  
  
  double p = c;
  
  for(k=0;k<KK;k++){
    
    tmp1 = r1*r2*1.0;
    
    tmp2 = as<double>( gsl_sf_poch(p-0.5,k));
    tmp3 = as<double>(gsl_sf_poch(p,2*k));
    tmp4 = as<double>(gsl_sf_fact(k));
    tmp5 = tmp2*tmp3*tmp4;		
    num1= as<double> (gsl_sf_pow_int(tmp1,k));
    num2=  as<double> (gsl_sf_hyperg_0F1(c+2*k,r1+r2));
    num3=  as<double> (gsl_sf_hyperg_0F1(c+2*k+1,r1+r2));
    //term[k]=(num1*num2*1.0/tmp5);
    temp=(num1*num2*1.0/tmp5);
    
    //tmp5 = tmp2*tmp3*tmp4;		
    
    temp=num1*(   ((2*k/d1)*num2) + ((d1/(2*c+4*k))*num3   )   )*1.0/tmp5;
    S +=temp ;
    if(  temp< MaxTolarableError) {
      /* terminate the loop using break statement */
      break;
    }
  }
  
  return S;
}



// [[Rcpp::export]]
double partial_d1_partial_d1_hyper_2by2_R(double c, NumericVector D, double MaxTolarableError=0.0000000000001){
  int KK = 30;
  double S=0.0;
  double temp,tmp1,tmp2,tmp3,tmp4,tmp5,num1,num2, num3,num4;
  int k=0;
  
  double r1=D[0];
  double r2=D[1];
  
  Rcpp::Environment gsl_env("package:gsl");
  Rcpp::Function gsl_sf_hyperg_0F1 = gsl_env["hyperg_0F1"];
  Rcpp::Function gsl_sf_fact= gsl_env["fact"];
  Rcpp::Function gsl_sf_poch= gsl_env["poch"];
  Rcpp::Function gsl_sf_pow_int= gsl_env["pow_int"];
  
  
  double d1 = 2*sqrt(r1);
  //double d2 = 2*sqrt(r2);
  
  
  double p = c;
  
  for(k=0;k<KK;k++){
    
    tmp1 = r1*r2*1.0;
    
    tmp2 = as<double>( gsl_sf_poch(p-0.5,k));
    tmp3 = as<double>(gsl_sf_poch(p,2*k));
    tmp4 = as<double>(gsl_sf_fact(k));
    tmp5 = tmp2*tmp3*tmp4;		
    num1= as<double> (gsl_sf_pow_int(tmp1,k));
    num2=  as<double> (gsl_sf_hyperg_0F1(c+2*k,r1+r2));
    num3=  as<double> (gsl_sf_hyperg_0F1(c+2*k+1,r1+r2));
    num4=  as<double> (gsl_sf_hyperg_0F1(c+2*k+2,r1+r2));
    //term[k]=(num1*num2*1.0/tmp5);
    
    temp=num1*(  ((2*k)*(2*k-1)/(d1*d1)* num2    ) + ((4*k+1)/(2*p+4*k) *num3) +  ((d1*d1)/((2*p+4*k)*(2*p+4*k+2))*num4)  )*1.0/tmp5;
    S += temp;
    if(  temp< MaxTolarableError) {
      /* terminate the loop using break statement */
      break;
    }
    
    
  }
  
  return S;
}





// [[Rcpp::export]]
double partial_d2_hyper_2by2_R(double c, NumericVector D, double MaxTolarableError=0.0000000000001){
  int KK = 30;
  double S=0.0;
  double temp,tmp1,tmp2,tmp3,tmp4,tmp5,num1,num2, num3;
  int k=0;
  
  double r1=D[0];
  double r2=D[1];
  
  Rcpp::Environment gsl_env("package:gsl");
  Rcpp::Function gsl_sf_hyperg_0F1 = gsl_env["hyperg_0F1"];
  Rcpp::Function gsl_sf_fact= gsl_env["fact"];
  Rcpp::Function gsl_sf_poch= gsl_env["poch"];
  Rcpp::Function gsl_sf_pow_int= gsl_env["pow_int"];
  
  //double d1 = 2*sqrt(r1);
  double d2 = 2*sqrt(r2);
  
  
  double p = c;
  
  for(k=0;k<KK;k++){
    
    tmp1 = r1*r2*1.0;
    
    tmp2 = as<double>( gsl_sf_poch(p-0.5,k));
    tmp3 = as<double>(gsl_sf_poch(p,2*k));
    tmp4 = as<double>(gsl_sf_fact(k));
    tmp5 = tmp2*tmp3*tmp4;		
    num1= as<double> (gsl_sf_pow_int(tmp1,k));
    num2=  as<double> (gsl_sf_hyperg_0F1(c+2*k,r1+r2));
    num3=  as<double> (gsl_sf_hyperg_0F1(c+2*k+1,r1+r2));
    // num4=  as<double> (gsl_sf_hyperg_0F1(c+2*k+2,r1+r2));
    //term[k]=(num1*num2*1.0/tmp5);
    
    
    temp= num1*(   ((2*k/d2)*num2) + ((d2/(2*p+4*k))*num3)   )*1.0/tmp5;
    S += temp;
    if(  temp< MaxTolarableError) {
      /* terminate the loop using break statement */
      break;
    }
    
  }
  
  return S;
}




// [[Rcpp::export]]
double partial_d2_partial_d2_hyper_2by2_R(double c, NumericVector D, double MaxTolarableError=0.0000000000001){
  
  int KK = 30;
  double S=0.0;
  double temp,tmp1,tmp2,tmp3,tmp4,tmp5,num1,num2, num3,num4;
  int k=0;
  
  double r1=D[0];
  double r2=D[1];
  
  Rcpp::Environment gsl_env("package:gsl");
  Rcpp::Function gsl_sf_hyperg_0F1 = gsl_env["hyperg_0F1"];
  Rcpp::Function gsl_sf_fact= gsl_env["fact"];
  Rcpp::Function gsl_sf_poch= gsl_env["poch"];
  Rcpp::Function gsl_sf_pow_int= gsl_env["pow_int"];
  
  
  //double d1 = 2*sqrt(r1);
  double d2 = 2*sqrt(r2);
  
  
  double p = c;
  
  for(k=0;k<KK;k++){
    
    tmp1 = r1*r2*1.0;
    
    tmp2 = as<double>( gsl_sf_poch(p-0.5,k));
    tmp3 = as<double>(gsl_sf_poch(p,2*k));
    tmp4 = as<double>(gsl_sf_fact(k));
    tmp5 = tmp2*tmp3*tmp4;		
    num1= as<double> (gsl_sf_pow_int(tmp1,k));
    num2=  as<double> (gsl_sf_hyperg_0F1(c+2*k,r1+r2));
    num3=  as<double> (gsl_sf_hyperg_0F1(c+2*k+1,r1+r2));
    num4=  as<double> (gsl_sf_hyperg_0F1(c+2*k+2,r1+r2));
    
    
    temp= num1*(  ((2*k)*(2*k-1)/(d2*d2)*num2) + ((4*k+1)/(2*p+4*k) *num3) +  ((d2*d2)/((2*p+4*k)*(2*p+4*k+2))*num4)  )*1.0/tmp5;
    S += temp;
    if(  temp< MaxTolarableError) {
      /* terminate the loop using break statement */
      break;
    }
    
  }
  
  return S;
}




// [[Rcpp::export]]
double  partial_d1_partial_d2_hyper_2by2_R(double c, NumericVector D, double MaxTolarableError=0.0000000000001){
  
  
  int KK = 30;
  double S=0.0;
  double temp,tmp1,tmp2,tmp3,tmp4,tmp5,num1,num2, num3,num4;
  int k=0;
  
  double r1=D[0];
  double r2=D[1];
  
  Rcpp::Environment gsl_env("package:gsl");
  Rcpp::Function gsl_sf_hyperg_0F1 = gsl_env["hyperg_0F1"];
  Rcpp::Function gsl_sf_fact= gsl_env["fact"];
  Rcpp::Function gsl_sf_poch= gsl_env["poch"];
  Rcpp::Function gsl_sf_pow_int= gsl_env["pow_int"];
  
  double d1 = 2*sqrt(r1);
  double d2 = 2*sqrt(r2);
  
  
  double p = c;
  
  for(k=0;k<KK;k++){
    tmp1 = r1*r2*1.0;
    
    tmp2 = as<double>( gsl_sf_poch(p-0.5,k));
    tmp3 = as<double>(gsl_sf_poch(p,2*k));
    tmp4 = as<double>(gsl_sf_fact(k));
    tmp5 = tmp2*tmp3*tmp4;		
    num1= as<double> (gsl_sf_pow_int(tmp1,k));
    num2=  as<double> (gsl_sf_hyperg_0F1(c+2*k,r1+r2));
    num3=  as<double> (gsl_sf_hyperg_0F1(c+2*k+1,r1+r2));
    num4=  as<double> (gsl_sf_hyperg_0F1(c+2*k+2,r1+r2));
    
    
    temp= num1*(  ((2*k)*(2*k)/(d1*d2)*num2) +( (k*(d1/d2 + d2/d1))/(p+2*k) *num3) +  ((d1*d2)/((2*p+4*k)*(2*p+4*k+2))*num4)  )*1.0/tmp5;
    S += temp;
    if(  temp< MaxTolarableError) {
      /* terminate the loop using break statement */
      break;
    }
  }
  
  return  S;
}

//}

