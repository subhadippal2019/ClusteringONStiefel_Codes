sliceDPM_Normal<-function(y, idel,maxiter,inc,lambda,gprior,s,c){

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
  n = length(y)
  
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
  mu = c(mu1,mu2,mu3)
 
  # number of classes
  nc = inc; 
  
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
    ## sample mu
    for (j in 1:nc)
    {
      ind = which(delta==j)
      xij = sum(y[delta[ind]])
      mj = length(ind)
      nmean = xij*lambda/(mj*lambda+s)
      nvar = 1/(mj*lambda+s)
      mu[j] = rnorm(1,nmean,sqrt(nvar))
    }
  #mu[1:3]=c(mu1,mu2,mu3)
   #mu=c(mu1,mu2,mu3)
    #
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
        mu[nc+1] = rnorm(1,0,sqrt(1/s))
        nc = nc+1
        mv = 1-v
      }
    }
   #nc=3
    ## sample delta
    for (i in 1:n)
    {
      probs = c()
      lambdapi = sqrt(lambda/2/pi)
      sj=c()
      kk=0
      for (j in (1:nc))
      {
        if (w[j]>u[i])
        { kk=kk+1
          probs[kk] = log(lambdapi)+ (-lambda*(y[i]-mu[j])^2/2)
          sj=c(sj,j)
        }
      }
      totprob = log(sum(exp(probs-max(probs)))) +max(probs)
      
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
  
  ret_list <- list(mu,w,delta,nc)
  return(ret_list) 
}
  