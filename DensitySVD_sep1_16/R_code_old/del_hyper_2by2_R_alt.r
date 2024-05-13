del_hyper_2by2_R_alt <- function(p,D){
  
  dyn.load("./C_code/hyper_2by2_R.so")
  
  eigenValues=D^2/4;
  
  d1 = D[1]
  d2 = D[2]
  
  KK = 30
  
  S = 0.0
  
  for(i in 1:KK){
    d = 0.0
    a = 1.5
    
    
    
    hyper0F1_val_a = .C("hyper_2by2_R",a,eigenValues,d)[[3]]
  
    d = 0.0
    a = 1.5+1
    hyper0F1_val_a_plus_one = .C("hyper_2by2_R",a,eigenValues,d)[[3]]
  
    val = (2.0/d1)*hyper0F1_val_a + (d1/2.0)*hyper0F1_val_a_plus_one
  }
  
  return (S)
}