#maxValueFind_for_d1<-function(d2,S11,N){
  
#  res = nlm(f=log_f_d1,d2,d2,S11,N)
  
#  res$minimum
  
#}

log_f_d1 <- function(d1,d2,S11,N){
  
  #dyn.load("./C_code/hyper_2by2_R.so")
  
  D=c(d1,d2)
  eigenValues = D^2/4
  d=0.0
  
  hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
  
  val = (d1*S11 - N*log(hyper0F1_val))
  
  return(val)
  
}