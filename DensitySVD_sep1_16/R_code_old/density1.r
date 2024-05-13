density1<-function(y,D,S,dimIndex,N,bLog){

  ## y=(1:4000)/1000
  ##out=density1(y,D=c(7,5.4),S=c(14,15),1,20,1)
  ##plot(y,exp(out-max(out)),type='l')
  #dyn.load("./C_code/hyper_2by2_R.so")
  dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  val=0*seq(1:length(y));
  
  for (i in  1:length(y)){
    
  x=y[i];
  alpha=1;
  beta=0;
  D[dimIndex]=x;
  eigenValues=D^2/4;
  d = 0.0
  if(bLog==1){
    hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
    val[i]=(alpha-1)*log(x)+(x*(S[dimIndex]-beta))- N * log(hyper0F1_val)
  }else{
    hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
    val[i]=x^(alpha-1)*exp(x*(S[dimIndex]-beta))/(hyper0F1_val)^N
  } 
  
  }
  #write.table(val,file="./density1.txt")
  return(val)  
  
}