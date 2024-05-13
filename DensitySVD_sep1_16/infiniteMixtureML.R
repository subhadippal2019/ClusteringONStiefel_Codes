update_z <- function(d,SB_alpha,SB_beta,z){
  #nc = length(unique(d))
  nc = length(z)
  a = b = z = rep(0,nc)
  for(j in 1:nc){
    a[j] = SB_alpha + sum(d==j)
    b[j] = SB_beta + sum(d>j)
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

add_new_cluster <-function(nc,SB_alpha,SB_beta,z,curr_param,hyper){
  
  z[nc] = rbeta(1,SB_alpha,SB_beta) 
  w = compute_w_from_z(z)
  #cat("z = ",z,"\n")
  
  M_param = hyper$G 
  #print(M_param)
  M = rmf.matrix(M_param) 
  ### TBD: make the max absolute value for columnwise vector positive 
  
  #### this should automatically make V with proper sign
  M = change_sign_M(M)
  
  #### conditional of V|M,D for Gibbs
  V_param = hyper$H 
  V = rmf.matrix(V_param)
  
  D = diag(sort(rgamma(2,hyper$alpha,hyper$beta),decreasing = T))
  
  curr_param$M = abind(curr_param$M,M)
  curr_param$D = abind(curr_param$D,D)
  curr_param$V = abind(curr_param$V,V)
  curr_param$F = abind(curr_param$F,M%*%D%*%V)
  cat("nc = ",nc,"\n")
  
  res = NULL
  res$z = z
  res$w = w
  res$curr_param = curr_param
  return(res)
}


### nonparametric mixture modeling of Matrix Langevin with fixed clusters
infiniteMixtureML <- function(data,nc,true_nc=-1,max_iter=2,run_id=1,output_file){
  
  ### j cluster index
  ### i data index
  ### y - 1xn vector, n is the number of data points
  ### n total number of data
  ### inc - the initial number of classes
  ### s - scaler for estimating mu, mu's prior is N(0,1/s)
  ### M - scaler hyperparameter for DP
  
  # number of classes
  
  ###### for stick-breaking ######
  SB_alpha = 0.5 #### to increase number of clusters decrease SB_alpha
  M = 1
  SB_beta = M
  ################################
  
  
  #library(gtools)
  #library(plyr)
  set.seed(43185)
  print(date())
  init_run()
  
  N = dim(data)[[3]]
  
  MCMC_sample = vector("list", max_iter)
  
  data_with_init_with_MCMC_samples = NULL	
  data_with_init_with_MCMC_samples$data = data 
  #MCMC_output_file = sprintf("MCMC_sample_%s_%s_%s.RData",nc,run_id,true_nc)
  #init_param_output_file = sprintf("init_param_MLE_%s_%s_%s.RData",nc,run_id,true_nc)
  ########### hyper parameters #######
  ### need to select empirical prior
  #   hyper=NULL
  #   hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
  #   hyper$H = (matrix(c(0,0,0,0),ncol = 2))
  #   hyper$alpha = 1.0
  #   hyper$beta = 0.0
  #   hyper$dir_alpha = rep(0.0,nc)
  #   for (j in 1:nc)
  #     hyper$dir_alpha[j] = 1.0
  #   
  ####################################
  
  #### initialization ####
  #Rprof(interval = 0.02)
  ###curr_param = init_param_true(nc)
  init_param = init_param_from_MLE(data,nc)
  data_with_init_with_MCMC_samples$init_param = init_param
  #save(init_param_MLE=init_param,file=init_param_output_file)
  ########################
  print(date())
  ########################
  vague_prior = 1
  hyper = hyper_selection(nc,vague_prior,init_param)
  ########################
  
  curr_param = init_param
  cluster_assign_vec = curr_param$id_arr
 
  #d_tmp_val = matrix(rep(0,N*nc),ncol=nc)
  
  #### from slice #####
  # sample z from beta(1,M)  
  z = rbeta(nc,SB_alpha,SB_beta)
  w = compute_w_from_z(z)
  ########################
  
  
  ## Main loop
  for (iter in 1:max_iter) {
    
    ### sample phi
    #mu = update_phi(y,d,mu)
    curr_param = cluster_param_update(data, curr_param, cluster_assign_vec, cluster_id_vec=1:nc,hyper)
    
    ### sample z
    z = update_z(cluster_assign_vec,SB_alpha,SB_beta,z) ### see Muller et al book
    w = compute_w_from_z(z)
 
    ## sample u
    u = runif(N)
    for (i in 1:N){
      u[i] = u[i]*w[cluster_assign_vec[i]]
    }
    

    res = NULL
    
    ustar = min(u)
    while (sum(w[1:nc]) <= (1-ustar)){
      nc = nc+1
      #alpha = 0.1 to make the new cluster proportion small
      res = add_new_cluster(nc,SB_alpha,SB_beta,z,curr_param,hyper) #### s must be structure proposal distributio for M,D,V
      cat("called adding new cluster \n")
      z = res$z
      w = res$w
      curr_param = res$curr_param
    }
    
    #cat(length(w),"\t",length(mu),"\n")
    ## sample d
    #browser()
    for (i in 1:N){
      w_greater_u_i = which(w > u[i])
      L = length(w_greater_u_i)
      pr_vec = rep(0,L)
      for(l in 1:L){
        sel_cluster_id = w_greater_u_i[l]
        pr_vec[l] = w[sel_cluster_id]*dMatrixLangevin(data[,,i],curr_param$M[,,sel_cluster_id],curr_param$D[,,sel_cluster_id],curr_param$V[,,sel_cluster_id]) ### variance is currently 1 + Normality assumption
      }
      #cat(pr_vec,"\n")
      if(length(pr_vec) == 0){
        cat("Error!!!!\n")
      }
      if(length(pr_vec) > 1){
        cluster_assign_vec[i] = sample(w_greater_u_i,1,pr_vec,replace=FALSE)
      }else{
        cluster_assign_vec[i] = w_greater_u_i[1]
      }
    }
    
    ######################
    if(iter%%10 == 0){
      print(paste0("MCMC main iter = ",iter))
    }
    ##### gibbs step #####
    #     t1=Sys.time()
    #     clust_assign_prob_mat_0 = matrix(rep(0,N*nc),nrow=N)
    #     for (i in 1:N){
    #       clust_assign_prob_i = rep(0,nc)
    #       for (cluster_id in 1:nc){
    #         d_tmp_val[i,cluster_id] = dMatrixLangevin(data[,,i],curr_param$M[,,cluster_id],curr_param$D[,,cluster_id],curr_param$V[,,cluster_id])
    #         
    #         clust_assign_prob_i[cluster_id] = pi_vec[cluster_id]*d_tmp_val[i,cluster_id]
    #       }
    #       clust_assign_prob_mat_0[i,] = clust_assign_prob_i
    #       #clust_assign_prob_i = clust_assign_prob_i/sum(clust_assign_prob_i)
    #       Z[i] = sample(1:nc,1,FALSE,clust_assign_prob_i)
    #       #print(clust_assign_prob_i)
    #     }
    #     Sys.time()-t1
    #     
    ### compact way to run the above commented code
    #t2=Sys.time()
    #clust_assign_prob_mat = matrix(rep(0,N*nc),nrow=N)
    #for (cluster_id in 1:nc){
    #  d_tmp_val[,cluster_id] = apply(data,3,function(x) dMatrixLangevin(x,curr_param$M[,,cluster_id],curr_param$D[,,cluster_id],curr_param$V[,,cluster_id]))
    #}    
    #clust_assign_prob_mat = d_tmp_val*t(matrix(rep(pi_vec,N),ncol=N))
    #Z = apply(clust_assign_prob_mat,1,function(x) sample(1:nc,1,FALSE,x))
    #Sys.time()-t2
    
    
    #cluster_assign_vec = Z
    #write.table(d_tmp_val,file="d_tmp_val.txt")
    #tmp_table = as.matrix(table(Z))
    #cluster_id_cnt = cbind(as.integer(row.names(tmp_table)),tmp_table)
    #cluster_id_cnt = rep(0,length(hyper$dir_alpha))
    
    #for (cluster_id in 1:nc){
    #  cluster_id_cnt[cluster_id] = sum(cluster_assign_vec == cluster_id)
    #}
    ## as an alternative to previous 3 statements
    cluster_id_cnt = sapply(1:nc, function(x) sum(cluster_assign_vec == x))
    
    #print(hyper$dir_alpha)
    print(cluster_id_cnt)
    #pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)
    
    #### for all the clusters ####
    #if(iter == 1){
    #  load("ML_dataset.Rdata")
    #  curr_param = L$curr_param
    #  cluster_assign_vec = L$clust
    #}
    #curr_param = cluster_param_update(data, curr_param, cluster_assign_vec, cluster_id_vec=1:nc,hyper) 
    
    ### need to update with updated cluster_assign_vec
    ##################################
    
    print(paste0("Mixture modeling: MCMC iteration update done ",iter))
    MCMC_sample[[iter]]$curr_param = curr_param
    MCMC_sample[[iter]]$pi_vec = NULL
    
    if(iter%%10==0){
      data_with_init_with_MCMC_samples$MCMC_sample = MCMC_sample 
      save(data_with_init_with_MCMC_samples, file = output_file)
    }
  }#### max_iter
  
  
  
  #G = NULL
  #G$curr_param = curr_param
  #G$clust = cluster_assign_vec
  #print(G$curr_param$D[,,1])
  
  #return(G)
  #save(D1_1,D2_1,D1_2,D2_2,D1_3,D2_3,file="D1_D2.Rdata")
  data_with_init_with_MCMC_samples$MCMC_sample = MCMC_sample 
  save(data_with_init_with_MCMC_samples, file = output_file)
  #save("MCMC_sample", file = MCMC_output_file)
  #Rprof(NULL)
  #out=summaryRprof()
  
  
  return(data_with_init_with_MCMC_samples)
  #return(out)
}


dMatrixLangevin <- function(X,M,D,V){
  
  #dyn.load("./C_code/functions_hyper_2by2_R.so")
  
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

init_run <-function(){
  source("utility.R")
  load_src_libs()
  
}

hyper_selection <-function(nc,vague_prior,init_param){
  
  hyper=NULL
  
  if(vague_prior == 1){
    print(paste0("using vague hyper-parameters"))
    
    hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
    hyper$H = (matrix(c(0,0,0,0),ncol = 2))
    hyper$alpha = 1.0
    hyper$beta = 0.15 ###### note that, min 5 data made a cluster 5*0.01 = 0.05, 
    ###### for N==1 if S>0.99 then distribution of d1 hv a heavy tail 
    
    hyper$dir_alpha = rep(0.0,nc)
    for (j in 1:nc)
      hyper$dir_alpha[j] = 1.0
  }else{ ## emperical
    
    print(paste0("using empirical hyper-parameters"))
    
    hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
    hyper$H = (matrix(c(0,0,0,0),ncol = 2))
    hyper$alpha = 1.0
    
    mean_D1_D2 = sum(init_param$D)/(2*nc)
    beta_mean = 1/mean_D1_D2
    
    hyper$beta = beta_mean
    
    hyper$dir_alpha = rep(0.0,nc)
    for (j in 1:nc)
      hyper$dir_alpha[j] = 1.0
  }
  
  
  hyper$debug = 0
  
  return(hyper)
  
}
