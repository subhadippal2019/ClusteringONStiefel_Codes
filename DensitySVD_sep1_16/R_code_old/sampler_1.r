library(rstiefel)
source("rdensity1.r")

N = 100

F = (matrix(c(1,2,3,4,5,6),ncol = 2))
sum_X = matrix(c(0,0,0,0,0,0),ncol=2)

for (i in 1:N)
{
  sum_X = sum_X + rmf.matrix(F) 
}
####################################
#### inference
####################################

G = (matrix(c(0,0,0,0,0,0),ncol = 2))
H = (matrix(c(0,0,0,0),ncol = 2))

out = svd(F)

D = out$d
M = out$u
V = out$v

#S = t(M)%*%sum_X%*%V

#nBin=100
#sampleSize=1
#dimIndex=1
#plot(exp(out-max(out)))

max_iter = 200

for (iter in 1:max_iter)
{
  #### conditional of D for Gibbs
  S = t(M)%*%sum_X%*%V
  ### for D[1] = d_1
  d1 = rdensity1(D,S,dimIndex=1,N,nBin=100,sampleSize=1)
  D[1] = d1
  ### for D[2] = d_2
  d2 = rdensity1(D,S,dimIndex=2,N,nBin=100,sampleSize=1)
  D[2] = d2
  
  #### conditional of M for Gibbs
  M_param = G + sum_X%*%t(V)%*%diag(D)
  M = rmf.matrix(M_param) 
  ### TBD: make the max absolute value for columnwise vector positive 
  
  #### conditional of V for Gibbs
  V_param = H + diag(D)%*%t(M)%*%sum_X
  V = rmf.matrix(V_param) 
  ### TBD: make the max absolute value for columnwise vector positive 
  
  
}


