library(rstiefel)
M_hyper_arr = array(c(1,0,0,0,1,0,1,0,0,0,-1,0,0,-1,0,0,0,-1), c(3, 2, 30))


M = array(rep(0,3*2*30), c(3, 2, 30))
set.seed(1)
nClust=1
for(ik in 1:nClust){
  M[,,ik] = rmf.matrix(M_hyper_arr[,,ik]) 
}

D0= diag(c(8,6),2,2)
D0= diag(c(10,7),2,2)
  #change D0 to get different simulation settings
#D = array(rep(0,2*2*30), c(2, 2, 30))
#for(ik in 1:nClust){
 # D[,,ik] = matrix(c(5+3*ik,0,0,3+3*ik))
#}

V = array(rep(0,2*2*30), c(2, 2, 30))
for(ik in 1:nClust){
  V[,,ik] = matrix(c(1,0,0,1))
}
# So V does not vary among clusters? - R.
##########################
N=2000
data = array(rep(0,3*2*N), c(3, 2, N)); 
## 0.2 0.4 0.4

#generating data-R
dl = rep(0,N)  # true cluster
for(i in 1:N){
  r = rmultinom(1,1,c(1/3,1/3,1/3))
  id = 1 # which(r==1)
  F = M[,,id]%*%D0%*%V[,,id]
  
  data[,,i] = rmf.matrix(F) 
  dl[i] = id
}
y=data


dyn.load("./C_code/functions_hyper_2by2_R.so")
dm<- function(X,M,D,V){
  F = M %*% D %*% V
  D = diag(D)
  eigenValues=D^2/4;
  dRet = 0.0
  ### normalizing constant
  hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,dRet)[[3]]
  Y1 = .C("partial_d1_hyper_2by2_R",1.5,eigenValues,dRet)[[3]]
  #YY1 = .C("partial_partial_hyper_2by2_R",1.5,eigenValues,dRet)[[3]]
  
  ### compute density value
nn=dim(X)[3]
ml=0; Sm=0;
for(ii in 1:nn){
  Sm=Sm+X[,,ii]
  #ml=ml+ (sum(diag(t(F)%*%X[,,ii]))) -log( hyper0F1_val )
  # do ml for just log likelihood.
}
S = t(M)%*%Sm%*%V
ly1 = S[1] - nn*Y1/hyper0F1_val
  #ML_density =exp(sum(diag(t(F)%*%X))) / hyper0F1_val
  return(ly1)
}


ik=1
M1= M[,,ik]
V1=V[,,ik]


#Compute partial log derivative for a range of d1 values
lg<-function(y,M1, D0, V1,es=.01){
ll=c()
kl=1
d2=D0[2,2]
kr=seq(d2,35,by=es)
Dx=D0
for  (d1 in kr ){
  Dx[1,1]=d1
  ll[kl]= dm(y, M1,Dx,V1)
  kl=kl+1
}
if(max(ll)<0) (ans= d2)
else{
  # go for mode finding
lgl=cbind(ll,kr)
  ll=lgl[,1]
  kr=lgl[,2]
#First step: Get two values with unequal signs
dx0=sample(1:length(kr),1)
dx1=dx0
sg0=sign(ll[dx0])
sgn=sg0
while  (sgn == sg0){
  dx0=dx1
  sg0=sign(ll[dx0])
  if (sgn== 1) dx1= sample(dx0:length(kr),1)
  if(sgn== -1) dx1=sample(1:dx0,1)
  sgn=sign(ll[dx1])
  }
#Second step: Use bisection method
dx0=min(dx0,dx1)
dx1=max(dx0,dx1)
#dr1=ll[1:(length(kr)-1)]
#dr2=ll[2:length(kr)]
#er=abs(dr1-dr2)
eps=.001
 while  ((abs(dx0-dx1)>1)) {
dx=ceiling((dx0+dx1)/2)   
if(sign(ll[dx]) >0)  dx0=dx
if(sign(ll[dx]) < 0)  dx1=dx
 }
ans=kr[dx0]
}
  return(ans)
}


