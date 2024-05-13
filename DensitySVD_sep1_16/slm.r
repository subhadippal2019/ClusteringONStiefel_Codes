source("dmode.r")
nClust = 3
M_hyper_arr = array(c(1,0,0,0,1,0,1,0,0,0,-1,0,0,-1,0,0,0,-1), c(3, 2, 30))
#30 arbitrary high number; to accomodate new clusters in sampling-R
# Since we are not sampling parameters, the 3 M, V, D s are fixed to simulation truth 
#throughout-R


#M_hyper_arr = array(c(1,0,0,0,1,0,1,0,0,0,-1,0,1/sqrt(2),0,1/sqrt(2),0,1,0), c(3, 2, 3))  
### mean matrix generation nClust matrices

# Generating parameters-R
M = array(rep(0,3*2*30), c(3, 2, 30))
set.seed(1)
for(i in 1:nClust){
  M[,,i] = rmf.matrix(M_hyper_arr[,,i]) 
}

D = array(rep(0,2*2*30), c(2, 2, 30))
for(i in 1:nClust){
  D[,,i] = matrix(c(5+3*i,0,0,3+3*i))
}

V = array(rep(0,2*2*30), c(2, 2, 30))
for(i in 1:nClust){
  V[,,i] = matrix(c(1,0,0,1))
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
  id = which(r==1)
  F = M[,,id]%*%D[,,id]%*%V[,,id]
  
  data[,,i] = rmf.matrix(F) 
 dl[i] = id
}
y=data


dyn.load("./C_code2/functions_hyper_2by2_R.so")
dmg<- function(X,M,D,V){
  F = M %*% D %*% V
  D = diag(D)
  eigenValues=D^2/4;
  dRet = 0.0
  ### normalizing constant
  hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,dRet)[[3]]
  
  ### compute density value
  ML_density =exp(sum(diag(t(F)%*%X))) / hyper0F1_val
  
  return(ML_density)
  
}


# Other hyperparms
epsilon = 0.5
s = 0.01
gprior = c(0.1,0.1)
# initial number of classes
inc = 3
# max number of iterations
maxiter = 600
# fix hyperparameter c
c = 0.5

# Main function-R
#sldpm<-function(y, idel,maxiter,inc,gprior,s,c){

  set.seed(123)
  #function [mu,w,delta,nc] = sliceDP(y,maxiter,inc,lambda,gprior,s,c)
  # slice sampler for the dirichlet mixture model
  # for one dimensional normal distribution with fixed sigma
  # Input:
  #     y - 1xn vector, n is the number of data points
  #     maxiter - the number of max iteration for Gibbs sampler
  #     inc - the initial number of classes
  #     lambda - sigma^-2 for normal distribution
  #     gprior - 1x2 vector the gamma prior for c
  #     s - scaler for estimating mu, mu's prior is N(0,1/s)
  #     c - scaler hyperparameter for DP
  # Output:
  #     mu - 1xm mean for each normal distribution
  #     w - 1xm weight for each component
  #     delta - 1xn class assignment for each data
  #     nc - final number of classes
  
  ## Initialization
  n = N # length(y)
  
  # sample v from beta(1,c)
  v = rbeta(inc,1,c)
  
  # compute w
  mv = 1-v
  w = v
  for (j in 2:length(v))
  {
    w[j] = v[j]*prod(mv[1:j-1])
  }
  
  
  
  # random initialize delta
  delta=ceiling(inc*runif(n))
  #delta = idel # ceiling(inc*runif(n))
  
  # initial mu

 
  # number of classes
  nc = inc; 
  
  
  hyper=NULL
  hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
  hyper$H = (matrix(c(0,0,0,0),ncol = 2))
  hyper$alpha = 1.0
  hyper$beta = 0.0
  hyper$dir_alpha = rep(0.0,nc)
  for (j in 1:nc)
    hyper$dir_alpha[j] = 1.0
  
  

  for(i in 1:nClust){
    M[,,i] = rmf.matrix(M_hyper_arr[,,i]) 
    ij=i+4
    D[,,i] = matrix(c(5+3*ij,0,0,3+3*ij))
  }
  
  
  
  ## Main loop
  for (iter in 1:maxiter) 
  {
    ## sample u
    u = runif(n)
    for (i in 1:n)
    {
      u[i] = u[i]*w[delta[i]]
    }
    # print(paste0("u ",u))
    ## sample theta
    #Suppress this loop  to fix theta
    for (ic in 1:nc){
      d1=lg(y[,,which(delta==ic)], M[,,ic], D[,,ic], V[,,ic])
      D[,,ic][1,1]=d1
      
      X = y[,,which(delta==ic)]
      nj =sum(delta==ic)
      sum_X = X[,,1]*0
      for (ij in 1:nj){
        sum_X = sum_X + X[,,ij]
      }
      
      M_param = hyper$G + sum_X%*%t(V[,,ic])%*%(D[,,ic])
      M[,,ic] = rmf.matrix(M_param) 
      ### TBD: make the max absolute value for columnwise vector positive 
      #### conditional of V|M,D for Gibbs
      V_param = hyper$H + (D[,,ic])%*%t(M[,,ic])%*%sum_X
      V[,,ic] = rmf.matrix(V_param) 
      ### TBD: make the max absolute value for columnwise vector positive 
      #### update the current parameter set
    }

    
    
    # print(paste0("mu ",mu))
    ## sample v and w
    # for j <= k*
    # update v
    # print(paste0("upto here nc = ",nc," dim v = ",length(v)))
    jk=sample(1:nc,1)
    for (j in c(jk))
    {
      # print(paste0("j ",j))
      ind = which(delta==j)
      # print(paste0("ind ",ind))
      if(length(ind) > 0)
      {
        # print(paste0("delta ",delta))
        alphaj = max(u[ind])
        mv = 1-v
        alphaj = alphaj/prod(mv[1:j-1])
        if(alphaj > 1)
        { #print(paste0("alphaj ",alphaj," j ",j," u ", u," v ",v,"\n"))
          }
      }else{
        alphaj = 0
      }
      # print(paste0("delta ", delta))
      ind = which(delta>j)
      # print(paste0("ind ",ind))
      if(length(ind) > 0)
      {
        ui = u[ind]
        # print(paste0("ui ",ui))
        vki = v[delta[ind]]
        mv[j] = 1
        prodki = rep(1,length(ind))
        for (l in 1:length(ind))
        {
          prodki[l] = prod(mv[1:delta[ind[l]]-1])
        }
        vki = vki*prodki
        betaj = 1-max(ui/vki)
      }else{
      #if(length(betaj) == 0)
      #{
        betaj = 1
      }
      
      zj = runif(1)
      v[j] = 1-((1-zj)*(1-alphaj)^c + zj*(1-betaj)^c)^(1/c)
      if(alphaj > 1)
       {
        # print(paste0("v ",v," alphaj ",alphaj," betaj ", betaj))
      }
      
    }
    # print(paste0("v ",v))
    # update w
    mv = 1-v
    w = v
    for (j in 2:length(v))
    {  
      w[j] = v[j]*prod(mv[1:j-1])
    }
    
    ustar = min(u)
    mustar = 1-ustar
  
    for (j in 1:nc)
    {
      # print(paste0("mustar ",mustar))
      # print(paste0("w ",w))
      if(sum(w[1:j])>mustar)
      {
        break
      }
    }
  
  
    if (j<nc)
    {  
      nc = j
    }else{
      while (sum(w[1:nc])<=mustar)
      {
        # sample more v from prior
        v[nc+1] = rbeta(1,1,c)
        w[nc+1] = v[nc+1]*prod(mv[1:nc])
        
  # These should be sampled from the base measure-R
        #whatever that is-R.
        gs=2
        M[,,nc+1] = rmf.matrix(M_hyper_arr[,,gs]) 
         # M[,,nc+1] = rmf.matrix(M_hyper_arr[,,nc+1]) 
        
   gh=runif(1,1,3)
   D[,,(nc+1)] = matrix(c(5+3*(gh),0,0,3+3*(gh)  ))
          #D[,,(nc+1)] = matrix(c(5+3*(nc+1),0,0,3+3*(nc+1)  ))
    
          V[,,(nc+1)] = matrix(c(1,0,0,1))
        
        
        #mu[nc+1] = rnorm(1,0,sqrt(1/s))
        nc = nc+1
        mv = 1-v
      }
    }
   nc=3
    ## sample delta
    for (i in 1:n)
    {
      probs = c()
      #lambdapi = sqrt(lambda/2/pi)
      sj=c()
      kk=0
      for (j in (1:nc))
      {
        if (w[j]>u[i])
        { kk=kk+1
          probs[kk] = log(dmg(y[,,i],M[,,j],D[,,j],V[,,j])   )         ##log(lambdapi)+ (-lambda*(y[i]-mu[j])^2/2)
          sj=c(sj,j)
        }
      }
      totprob = log(sum(exp(probs-max(probs)))) +max(probs)
      #log scale to prevent underflow-R
  #The following piece is alright, but I just used  sample() to make this shorter.-R 
      
#       r = totprob*runif(1)
#       maxp = probs[1]
#       cl = 1
#       while(r>maxp)
#       {
#         cl = cl+1
#         maxp = maxp+probs[cl]
#       }
    if(length(sj)>1) delta[i]=sample(sj,size=1,prob= exp(probs-totprob))              # delta[i] = cl
      if(length(sj)==1) delta[i]=sj
    }
    
    #mcmc ends
    
  }
  
 # ret_list <- list(w,delta,nc)
 # When the par-s will be sampled, it will also return cluster specific F, or M, D V 
 #  return(ret_list) 
#}

#a=sldpm(y, idel,maxiter,inc,gprior,s,c)
  