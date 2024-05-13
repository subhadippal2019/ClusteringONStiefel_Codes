## test file for slice sampling

## load data
### need to generate data
w1 = 0.3
w2 = 0.2
w3 = 0.5
mu1 = 1
mu2 = 0
mu3 = -1
sig = .1
lambda=sig^(-2)
#\sum w_i N(\mu_i,1)
nData=2000
set.seed(2)
Data = rep(0,nData)
dl=c()
for(i in 1:nData )
{
  r = runif(1)  
  if(r <= w1) {
    Data[i] = rnorm(1,mu1,sig) 
  dl[i]=1 }
  if(r > w1 && r <= (w1+w2)){
    Data[i] = rnorm(1,mu2,sig) 
  dl[i]=2 }
  if(r > (w1+w2)) {
    Data[i] = rnorm(1,mu3,sig) 
  dl[i]=3}
}

## set parameters
epsilon = 0.5
s = 0.01
gammaprior = c(0.1,0.1)

# initial number of classes
inc = 3
# max number of iterations
maxiter = 1000

# fix hyperparameter c
c = 0.5

source("sl.r")
# slice sampling the dirichlet mixture model
#dl is the cluster index vector.
A = sliceDPM_Normal(Data,dl,maxiter,inc,lambda,gammaprior,s,c)

