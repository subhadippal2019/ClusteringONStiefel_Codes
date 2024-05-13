dyn.load("functions_hyper_2by2_R.so")

V1 = NULL
V2 = NULL

D = NULL
K=30
X = seq(1,K,by=0.1)

a=1.5

D_hat_1 = 15
D_hat_2 = 6

#for (i in 1:length(X)){
  dRet = 0.0
  
  D[2] = D_hat_2
  D[1] = D_hat_1    
  
  eigenValues=D^2/4
    
  h = .C("hyper_2by2_R",cc=a,eigenValues,dRet)[[3]]
  h1 = .C("partial_d1_hyper_2by2_R",cc=a,eigenValues,dRet)[[3]]
  h2 = .C("partial_d2_hyper_2by2_R",cc=a,eigenValues,dRet)[[3]]
  
  h1/h
  h2/h
  
  V1[i] = h
  V2[i] = h1
  
#}

plot(X,V2/V1,type='l',col='red',xlim=c(0,K))
abline(0,1)
