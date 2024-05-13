rdensity1<-function(D,S,dimIndex,N,nBin=10,sampleSize=1){
  
  print(paste0("inside rdensity"))
  ##out=density1(y=(1:4000)/2000,D=c(0.7,5.4),S=c(14,15),1,10,1)
  ##plot(exp(out-max(out)))
  dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  #### y is the supprt of the distribution discretized in a grid
  #### if dimIndex = 1 then we have to sample for d_1 (> d_2)
  #### if dimIndex = 2 then we have to sample for d_2 (< d_1)
  if(length(dim(D)) == 2 ){
    D = diag(D)
  }
  if(length(dim(S)) == 2 ){
    S = diag(S)
  }
  if(nBin > 100)
  {
    nBin = 100
  }

  
  if(dimIndex == 1) #### d_1 should be in d_2 to \infty
  {
    d1_near_mode = NR_method_d1(D,S,dimIndex,N,max_iter=20) ##### replaced by NR method 
    diff_from_d2 = abs(D[2]-d1_near_mode) #### how far should we search; need a lower bound
    pseudoInfty = D[2] + 5*diff_from_d2   #### totally adhoc
    y = seq(D[2],pseudoInfty,1/nBin) 
  }
  if(dimIndex == 2) #### d_2 should be 0 to d_1
  {
    if(D[1] <= 1/nBin)
    {
      nBin = 100*nBin
    }
    y = seq(1/nBin,D[1],1/nBin) 
  }
  
  val=0*seq(1:length(y));
  
  for (i in  1:length(y)){
    
    x=y[i];
    alpha=1;
    beta=0;
    D[dimIndex]=x;
    eigenValues=D^2/4;
    d = 0.0
    bLog = 1 ### hard-coded to avoid problem
    if(bLog==1){
      hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
      val[i]=(alpha-1)*log(x)+(x*(S[dimIndex]-beta))- N * log(hyper0F1_val)
    }else{
      hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
      val[i]=x^(alpha-1)*exp(x*(S[dimIndex]-beta))/(hyper0F1_val)^N
    } 
    
  }
  #print(val)
  #write.table(val,file="./rdensity1.txt")
  prob=exp(val-max(val))
  plot(prob,type='l')
  
  S1=y[which(prob>0)]
  prob1=prob[which(prob>0)]
  #S1
  x=sample(S1, size=sampleSize, replace = TRUE, prob = prob1)
  
  return(x)  
  
}