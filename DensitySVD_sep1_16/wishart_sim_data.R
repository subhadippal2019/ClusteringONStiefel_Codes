#######
### try to generate data that could work with mixture of Matrix-Langevins 
### but may not work with mixture of wisharts
### 
#######
#### case 1: generate two sets of data with M_1,M_2 keeping D same

trace <-function(X){
  if(dim(X)[1] > 1)
    if(dim(X)[1] == dim(X)[2])
      return(sum(diag(X)))
}

dwishart <-function(x,n,V){
  p = dim(x)[2]
  log_density = (n-p-1)/2*log(det(x)) - trace(  solve(V)%*%x  )/2 - n*p/2*log(2) -
                n/2*log(det(V)) - p*(p-1)/4*log(pi) - sum(lgamma(n/2-(seq(1:p)-1)/2))     
  d = (log_density) 

  return(d)
}

drawEllipse <-function(M){

  if(sum(   (t(M)%*%M-diag(dim(t(M)%*%M)[2]))^2  )<.0001) {  A=M%*%diag(c(3,0.2))%*%t(M);
  }  else {A=M; print(A)}
  #if(!sum(t(M)%*%M-diag(dim(t(M)%*%M)[2]))^2) {  A=M}
  
  ctr    <- c(0, 0)                               # data centroid -> colMeans(dataMatrix)
  #A      <- matrix(c(2.2, 0.4, 0.4, 2.8), nrow=2) # covariance matrix -> cov(dataMatrix)
  #RR     <- chol(A)                               # Cholesky decomposition
  angles <- seq(0, 2*pi, length.out=200)          # angles for ellipse
  #ell    <- 1 * cbind(cos(angles), sin(angles)) %*% RR  # ellipse scaled with factor 1
  #ellCtr <- sweep(ell, 2, ctr, "+")               # center ellipse to the data centroid
  #plot(ellCtr, type="l", lwd=2, asp=1)            # plot ellipse
  #points(ctr[1], ctr[2], pch=4, lwd=2)            # plot data centroid

  eigVal  <- eigen(A)$values
  eigVec  <- eigen(A)$vectors
  eigScl  <- eigVec %*% diag(sqrt(eigVal))  # scale eigenvectors to length = square-root
  xMat    <- rbind(ctr[1] + eigScl[1, ], ctr[1] - eigScl[1, ])
  yMat    <- rbind(ctr[2] + eigScl[2, ], ctr[2] - eigScl[2, ])
  ellBase <- cbind(sqrt(eigVal[1])*cos(angles), sqrt(eigVal[2])*sin(angles)) # normal ellipse
  ellRot  <- eigVec %*% t(ellBase)                                          # rotated ellipse
  plot((ellRot+ctr)[1, ], (ellRot+ctr)[2, ], asp=1, type="l", lwd=2)
  matlines(xMat, yMat, lty=1, lwd=2, col="green")
  points(ctr[1], ctr[2], pch=4, col="red", lwd=3)
}

estWishartParam <-function(C1){

  N = dim(C1)[3]
  E_X_1 = apply(C1,c(1,2),'mean')
  V_X_1 = apply(C1,c(1,2),'var')
  V1 = diag(2)
  
  V1[1,1] = 0.5*V_X_1[1,1]/E_X_1[1,1]
  V1[2,2] = 0.5*V_X_1[2,2]/E_X_1[2,2]
  K1 = V1[1,1]*V1[2,2]
  R1 = V_X_1[1,2]/E_X_1[1,2]
  V1[1,2] = V1[2,1] = (R1+sqrt(R1*R1-4*K1))/2
  
  
  df = (E_X_1[1,1]+E_X_1[1,2]+E_X_1[2,2])/(V1[1,1]+V1[1,2]+V1[2,2])
  
  return(list(V=V1,df=df))
}

NR_method_for_wishart <-function(C1,n_old=2,MAX_ITER = 50){

  N = dim(C1)[3]
  S1 = mean(log(apply(C1, 3, 'det'))) # \frac{1}{N}\sum_{i=1}^{N}log(|X_i|)
  S2 = log(det(apply(C1, c(1,2), 'sum')))
  iter = 1
  
  while(iter < MAX_ITER){
    
    f_n_old = digamma(n_old/2) + digamma((n_old-1)/2) - 2* log(N*n_old/2) + S2 - S1
    f_prime_n_old = 0.5*(trigamma(n_old/2) + trigamma((n_old-1)/2)) - 2/n_old
  
    n_new = n_old - f_n_old/f_prime_n_old
    cat(n_new,"\n")
    n_old = n_new
    iter = iter+1
  }
  return(n_new)
}


N = 1000
C1 = array(rep(N*2*2),c(2,2,N))
C2 = array(rep(N*2*2),c(2,2,N))
C = array(rep(N*2*2),c(2,2,2*N))

Theta_1 = rep(0,N)
Theta_2 = rep(0,N)

for(i in 1:N){
  D_1 = (diag(c(6,2)+rep(runif(1,0,0.5),2)))/10
  theta = (pi/2) + runif(1,0,0.02)
  M_1 = matrix(c(sin(theta),cos(theta),-cos(theta),sin(theta)),nrow=2)
  C1[,,i] = M_1%*%D_1%*%t(M_1)  
  Theta_1[i] = theta
  
  D_2 = (diag(c(6,3)+rep(runif(1,0,10),2)))/10
  theta = (pi/4+ runif(1,0,0.02))
  M_2 = matrix(c(sin(theta),cos(theta),-cos(theta),sin(theta)),nrow=2)
  C2[,,i] = M_2%*%D_2%*%t(M_2)
  Theta_2[i] = theta
}

C[,,1:N] = C1
C[,,((N+1):(2*N))] = C2

par(mfrow=c(1,2))
drawEllipse(eigen(C1[,,200])$vectors)
drawEllipse(eigen(C2[,,200])$vectors)

par(mfrow=c(1,2))
drawEllipse((C1[,,200]))
drawEllipse((C2[,,200]))

n1 = NR_method_for_wishart(C1,3)
V1 = apply(C1,c(1,2),'mean')/n1

n2 = NR_method_for_wishart(C2,3)
V2 = apply(C2,c(1,2),'mean')/n2

library(MCMCpack)

dC1 = apply(C,3,function(x){dwish(x,n1,V1)})
dC2 = apply(C,3,function(x){dwish(x,n2,V2)})

dC3 = apply(C,3,function(x){dwishart(x,n1,V1)})
dC4 = apply(C,3,function(x){dwishart(x,n2,V2)})

dC1=dC3
dC2=dC4
plot(dC1,dC1,col='red',xlim = c(0,0.2), ylim=c(0,0.2))
par(new=T)
plot(dC2,dC2,col='blue',xlim = c(0,0.2), ylim = c(0,0.2))

hist(dC1[1:1000]/dC2[1:1000])
hist(dC1[1001:2000]/dC2[1001:2000])

plot(dC1[1:1000]/dC2[1:1000])
plot(dC1[1001:2000]/dC2[1001:2000])
 
which(dC1/dC2>1)
par(mfrow=c(1,2))
plot(dC1[1:1000]/dC2[1:1000]>1)
plot(dC1[1001:2000]/dC2[1001:2000]>1)
sum(dC1/dC2>1)

# simulation setup may me. same direction but different scale. Wishart will give different result. Non identification of fiber.

#### data generation from wishart to test how robust our estimators are
X = rWishart(10000,3.4,diag(c(0.2,1)))

X1 = X
X1 = apply(X,3,function(y){ z =  y+diag(rep(runif(1,0,0.5),2)); return (matrix(z,nrow=2))} )
X1 = array(X1,c(2,2,N))
n1 = NR_method_for_wishart(X1,3)
V1 = apply(X1,c(1,2),'mean')/n1

X2 = X
X2 = apply(X,3,function(y){ z =  y+diag(rep(runif(1,0,50),2)); return (matrix(z,nrow=2))} )
X2 = array(X2,c(2,2,N))
n2 = NR_method_for_wishart(X2,3)
V2 = apply(X2,c(1,2),'mean')/n2

