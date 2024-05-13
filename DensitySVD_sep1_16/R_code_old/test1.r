

require(rstiefel)
source("density1.r")

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

S = t(M)%*%sum_X%*%V

nBin = 100
MAX_X = 200
d1 = seq(D[2],MAX_X,1/nBin)

dimIndex = 1
bLog = 1

out1 = density1(y=d1,D,S,dimIndex,N,bLog)
out2 = density2(y=d1,D,S,dimIndex,N,bLog)

plot(d1,exp(out1-max(out1)))
par(new=TRUE)
plot(d1,exp(out2-max(out2)))
