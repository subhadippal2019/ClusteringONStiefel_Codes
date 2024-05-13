del_density1_alt<-function(y,D,S,dimIndex,N,bLog){
  
  ##out=density1(y=(1:6000)/200,D=c(0.7,5.4),S=c(14,15),1,10,1)
  ##plot(exp(out-max(out)))
  
  dyn.load("./C_code/del_hyper_2by2_R_alt.so")
  dyn.load("./C_code/hyper_2by2_R.so")
  
  val=0*seq(1:length(y));
  
  for (i in  1:length(y)){
    
    x=y[i];
    alpha=1;
    beta=0;
    D[dimIndex]=x;
    eigenValues=D^2/4;
    d = 0.0
    #     if(bLog==1){
    #       hyper0F1_val = .C("del_hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
    #       val[i]=(alpha-1)*log(x)+(x*(S[dimIndex]-beta))- N * log(hyper0F1_val)
    #     }else{
    #       hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
    #       val[i]=x^(alpha-1)*exp(x*(S[dimIndex]-beta))/(hyper0F1_val)^N
    #     } 
    #     
    Y1 = .C("del_hyper_2by2_R_alt",1.5,eigenValues,d)[[3]]
    hyper0F1_val = .C("hyper_2by2_R",1.5,eigenValues,d)[[3]]
    val[i] = S[dimIndex] - N*Y1/hyper0F1_val
  }
  #write.table(val,file="./density1.txt")
  return(val)  
  
}

del_density1_alt_1<-function(y,D,S,dimIndex,N,bLog){
  
  ##out=density1(y=(1:6000)/200,D=c(0.7,5.4),S=c(14,15),1,10,1)
  ##plot(exp(out-max(out)))
  
  #dyn.load("./C_code/del_hyper_2by2_R_alt.so")
  dyn.load("./C_code/hyper_2by2_R.so")
  
  val=0*seq(1:length(y));
  
  for (i in  1:length(y)){
    
    x=y[i];
    alpha=1;
    beta=0;
    D[dimIndex]=x;
    eigenValues=D^2/4;
    d = 0.0
    #     if(bLog==1){
    #       hyper0F1_val = .C("del_hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
    #       val[i]=(alpha-1)*log(x)+(x*(S[dimIndex]-beta))- N * log(hyper0F1_val)
    #     }else{
    #       hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
    #       val[i]=x^(alpha-1)*exp(x*(S[dimIndex]-beta))/(hyper0F1_val)^N
    #     } 
    #     
    #Y1 = .C("del_hyper_2by2_R_alt",1.5,eigenValues,d)[[3]]
    hyper0F1_val_1 = .C("hyper_2by2_R",1.5,eigenValues,d)[[3]]
    D[dimIndex]=x+0.001;
    eigenValues=D^2/4;
    hyper0F1_val_2 = .C("hyper_2by2_R",1.5,eigenValues,d)[[3]]
    
    Y1 = (hyper0F1_val_2-hyper0F1_val_1) / 0.001
    
    val[i] = S[dimIndex] - N*Y1/hyper0F1_val_1
  }
  #write.table(val,file="./density1.txt")
  return(val)  
  
}


partial_density1<-function(y,D,S,dimIndex,N,bLog){
  
  ##out=density1(y=(1:6000)/200,D=c(0.7,5.4),S=c(14,15),1,10,1)
  ##plot(exp(out-max(out)))
  
  dyn.load("./C_code/hyper_2by2_R.so")
  dyn.load("./C_code/partial_hyper_2by2_R.so")
  
  val=0*seq(1:length(y));
  
  for (i in  1:length(y)){
    
    x=y[i];
    alpha=1;
    beta=0;
    D[dimIndex]=x;
    eigenValues=D^2/4;
    d = 0.0
       
    Y1 = .C("partial_hyper_2by2_R",1.5,eigenValues,d)[[3]]
    hyper0F1_val = .C("hyper_2by2_R",1.5,eigenValues,d)[[3]]
    val[i] = S[dimIndex] - N*Y1/hyper0F1_val
  }
  #write.table(val,file="./density1.txt")
  return(val)  
  
}


partial_partial_density1<-function(y,D,S,dimIndex,N,bLog){
  
  ##out=density1(y=(1:6000)/200,D=c(0.7,5.4),S=c(14,15),1,10,1)
  ##plot(exp(out-max(out)))
  
  dyn.load("./C_code/hyper_2by2_R.so")
  dyn.load("./C_code/partial_hyper_2by2_R.so")
  
  val=0*seq(1:length(y));
  
  for (i in  1:length(y)){
    
    x=y[i];
    alpha=1;
    beta=0;
    D[dimIndex]=x;
    eigenValues=D^2/4;
    d = 0.0
    
    Y1 = .C("partial_hyper_2by2_R",1.5,eigenValues,d)[[3]]
    
    YY1 = .C("partial_partial_hyper_2by2_R",1.5,eigenValues,d)[[3]]
    
    hyper0F1_val = .C("hyper_2by2_R",1.5,eigenValues,d)[[3]]
    val[i] = - N*(YY1/hyper0F1_val - (Y1/hyper0F1_val)*(Y1/hyper0F1_val))
  }
  #write.table(val,file="./density1.txt")
  return(val)  
  #out2 = partial_partial_density1(y=(1:8000)/200,D,S,1,15,1)
}
