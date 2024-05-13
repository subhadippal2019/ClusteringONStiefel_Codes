
hyper_2by2_R_1 <-function(c,R){
  
  library(gsl)  
  
  ptm = proc.time()
  
  dyn.load("./C_code/functions_hyper_2by2_R.so")
  dRet = 0.0
  for(i in 1:100){
    hyper0F1_val = .C("hyper_2by2_R",cc=1.5,R,dRet)[[3]]
  }
  print(hyper0F1_val)
  print(proc.time() - ptm)
  
  ptm = proc.time()
  
  KK = 30
  S=0.0
  k=0
  
  r1=R[1]
  r2=R[2]
  

  for(i in 1:100){
  for(k in (0:(KK-1))){
    
    tmp1 = r1*r2*1.0;
    
    tmp2 = poch(c-0.5,k,give=FALSE,strict=TRUE)
    tmp3 = poch(c,2*k,give=FALSE,strict=TRUE)
    tmp4 = fact(k,give=FALSE,strict=TRUE)
    
    tmp5 = tmp2*tmp3*tmp4    
    
    S = S + pow_int(tmp1,k,give=FALSE,strict=TRUE)*hyperg_0F1(c+2*k,r1+r2,give=FALSE, strict=TRUE)*1.0/tmp5
    
  }  
  }
  print(proc.time() - ptm)
  print(S)
  return(S)
}   
