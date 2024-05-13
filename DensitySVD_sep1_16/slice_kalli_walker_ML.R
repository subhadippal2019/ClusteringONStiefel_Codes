

update_phi <- function(y, d, phi_old){
  nc = length(phi_old)
  phi = rep(0,nc)
  for(j in 1:nc){
    if(sum(d==j) > 0){
      phi[j] = mean(y[d==j])
    }else{
      phi[j] = phi_old[j]
    }
  }
  return(phi)
}

update_z <- function(d,alpha,beta,z){
  #nc = length(unique(d))
  nc = length(z)
  a = b = z = rep(0,nc)
  for(j in 1:nc){
    a[j] = alpha + sum(d==j)
    b[j] = beta + sum(d>j)
    z[j] = rbeta(1,a[j],b[j])
  }
  return(z)
}

compute_w_from_z <-function(z){
  # compute w
  mz = 1-z
  w = z
  if(length(z) > 1){
    for (j in 2:length(z)){
      w[j] = z[j]*prod(mz[1:(j-1)])
    }
  }else{
    cat("**** only 1 cluster left *****\n")
  }
  return(w)
}

add_new_cluster <-function(nc,alpha=1,beta,s,z,mu){
  z[nc] = rbeta(1,alpha,beta) 
  w = compute_w_from_z(z)
  #cat("z = ",z,"\n")
  mu[nc] = rnorm(1,0,1/s)
  res = NULL
  res$z = z
  res$w = w
  res$mu = mu
  return(res)
}

# Main function-R
#y, maxiter, inc, s, alpha, beta
slice_DPM_Normal <- function(y,maxiter=500,nc=1,s=0.1,alpha=0.5,beta=1,classification=NULL){
  
  res = NULL
  
  ## Initialization
  n = length(y)
  
  
  if(is.null(classification)){
    classification = ceiling(nc*runif(n))
  }
  ### j cluster index
  ### i data index
  ### y - 1xn vector, n is the number of data points
  ### n total number of data
  ### inc - the initial number of classes
  ### s - scaler for estimating mu, mu's prior is N(0,1/s)
  ### M - scaler hyperparameter for DP
  
  # number of classes
  
  # sample z from beta(1,M)
  z = rbeta(nc,alpha,beta)
  w = compute_w_from_z(z)
  
  # random initialize d_i
  d = classification
  #d = ceiling(nc*runif(n))

  # initial mu
  mu = rnorm(nc,0,1/s)  ### this is phi
  
  ## Main loop
  for (iter in 1:maxiter) {
    ### sample phi
    mu = update_phi(y,d,mu)
    
    ### sample z
    z = update_z(d,alpha,beta,z) ### see Muller et al book
    w = compute_w_from_z(z)
    
    ## sample u
    u = runif(n)
    for (i in 1:n){
      u[i] = u[i]*w[d[i]]
    }
  
    ## N = n
  
    ustar = min(u)
    while (sum(w[1:nc]) <= (1-ustar)){
      nc = nc+1
      #alpha = 0.1 to make the new cluster proportion small
      res = add_new_cluster(nc,alpha,beta,s,z,mu)
      cat("called adding new cluster \n")
      z = res$z
      w = res$w
      mu = res$mu
    }
    
    #cat(length(w),"\t",length(mu),"\n")
    ## sample d
    #browser()
    for (i in 1:n){
      w_greater_u_i = which(w > u[i])
      L = length(w_greater_u_i)
      pr_vec = rep(0,L)
      for(l in 1:L){
        k = w_greater_u_i[l]
        
        pr_vec[l] = w[k]*dnorm(y[i],mu[k],1) ### variance is currently 1 + Normality assumption
      }
      #cat(pr_vec,"\n")
      if(length(pr_vec) == 0){
        cat("Error!!!!\n")
      }
      if(length(pr_vec) > 1){
        d[i] = sample(w_greater_u_i,1,pr_vec,replace=FALSE)
      }else{
        d[i] = w_greater_u_i[1]
      }
    }
    
  }
  res$mu = mu
  res$w = w
  res$d = d

  return(res)
}
    

source("cluster_param_update.R")
#############################
set.seed(761583)
#############################
## data generation
N = 1000
pi_vec = c(0.2,0.2,0.2,0.2,0.2)
u = apply(rmultinom(N,1,pi_vec),2,function(x) which(x==1))

y = (u==1)*rnorm(N,3,1) + (u==2)*rnorm(N,-4,1) + (u==3)*rnorm(N,10,1) + (u==4)*rnorm(N,20,1) + (u==5)*rnorm(N,-10,1)

##############################
### initial clustering
library(mclust)
MC = Mclust(y,1:5)
cat("nc = ", dim(MC$z)[2],"\n")
#####################
### call slice sampler of DP


#id_tmp = (u==3) 
#unif_sel = rbinom(N,1,0.5)
#id_3 = which(id_tmp*unif_sel==1)
#id_4 = which(id_tmp*(1-unif_sel)==1)

#MC$classification[id_3] = 2
#MC$classification[id_4] = 3


maxiter = 500
inc = length(unique(MC$classification))
s = 0.03
alpha = 0.5
M = 1
beta = M

#res = slice_DPM_Normal(y, maxiter, inc, s, alpha,beta,MC$classification)
res = slice_DPM_Normal(y)

id = as.integer(names(table(res$d)))

cat(res$mu[id],"\n")
cat(res$w[id],"\n")
print(table(res$d))

######################