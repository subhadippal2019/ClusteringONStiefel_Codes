analyse_MCMC_sample <-function(MCMC_sample,L,perm_vec=c(1,2,3)){
  
  
  #load("ML_dataset.RData")
  
  nc = 3
  n_sam = length(MCMC_sample)
  D1_dist = matrix(rep(0,nc*n_sam),ncol=n_sam)
  D2_dist = matrix(rep(0,nc*n_sam),ncol=n_sam)
  
  V_dist = matrix(rep(0,nc*n_sam),ncol=n_sam)
  #M_dist = matrix(rep(0,nc*n_sam),ncol=n_sam)
  
  F_dist = matrix(rep(0,nc*n_sam),ncol=n_sam)
  
  #perm_vec = c(1,2,3)
  for (i in 1:n_sam){
    
    curr_param = MCMC_sample[[i]]$curr_param
    
    for(j in 1:nc){
      
      D1_dist[j,i] = curr_param$D[,,j][1,1] - L$curr_param$D[,,perm_vec[j]][1,1]
      D2_dist[j,i] = curr_param$D[,,j][2,2] - L$curr_param$D[,,perm_vec[j]][2,2]
    
      #V_dist[j,i] = 2 - sum(diag(t(curr_param$V[,,j])%*%L$curr_param$V[,,j]))
    
      #M_dist[j,i] = 2 - sum(diag(t(curr_param$M[,,j])%*%L$curr_param$M[,,j]))
      
      F_dist[j,i] = sum(abs(curr_param$M[,,j]%*%curr_param$D[,,j]%*%curr_param$V[,,j]-L$curr_param$M[,,perm_vec[j]]%*%L$curr_param$D[,,perm_vec[j]]%*%L$curr_param$V[,,perm_vec[j]]))
      V_dist[j,i] = (curr_param$V[,,j])[2,1]/(curr_param$V[,,j])[1,1]  
    } 
  }

  for(j in 1:nc){
    f_name = sprintf("clust_dist_%d.pdf",j)
    pdf(f_name)
    par(mfrow = c(2,2))
    plot(D1_dist[j,],type='l',col = 'red')
    plot(D2_dist[j,],type='l',col = 'red')
    plot(V_dist[j,],type='l')
    plot(F_dist[j,],type='l')
    dev.off()
  }
  

}




plot_mcmc_sample_D <-function(MCMC_sample,L){
  library(MASS)
  n_sam = length(MCMC_sample)
  D1=matrix(rep(0,0,2*n_sam),ncol=2)
  D2=matrix(rep(0,0,2*n_sam),ncol=2)
  D3=matrix(rep(0,0,2*n_sam),ncol=2)
  
  for(i in 1:n_sam){
    D1[i,] = diag(MCMC_sample[[i]]$curr_param$D[,,1])
    D2[i,] = diag(MCMC_sample[[i]]$curr_param$D[,,2])
    D3[i,] = diag(MCMC_sample[[i]]$curr_param$D[,,3])
  }
  
  pdf("cluster_1_D.pdf")
  f1 = kde2d(D1[,1], D1[,2], n = 200)
  persp(f1,theta = 60, phi = 30,border=NA,col = 'blue',axes=TRUE,box=TRUE,shade=0.5)
  contour(f1,col='blue')
  points(diag(L$curr_param$D[,,1])[1],diag(L$curr_param$D[,,1])[2],col = 'red',pch=16)
  dev.off()
  
  pdf("cluster_2_D.pdf")
  f1 = kde2d(D2[,1], D2[,2], n = 200)
  persp(f1,theta = 60, phi = 30,border=NA,col = 'blue',axes=TRUE,box=TRUE,shade=0.5)
  contour(f1,col='blue')
  points(diag(L$curr_param$D[,,2])[1],diag(L$curr_param$D[,,2])[2],col = 'red',pch=16)
  dev.off()
  
  pdf("cluster_3_D.pdf")
  f1 = kde2d(D3[,1], D3[,2], n = 200)
  persp(f1,theta = 60, phi = 30,border=NA,col = 'blue',axes=TRUE,box=TRUE,shade=0.5)
  contour(f1,col='blue')
  points(diag(L$curr_param$D[,,3])[1],diag(L$curr_param$D[,,3])[2],col = 'red',pch=16)
  dev.off()
  
}

cut_MCMC_samples <- function(MCMC_sample,start,end){
  
  cut_MCMC_sample = NULL
  j = 0
  for(i in start:end){
    j = j+1
    cut_MCMC_sample[[j]] = MCMC_sample[[i]] 
  }  
  return (cut_MCMC_sample)
}




calculate_DIC <-function(MCMC_sample,data,burn_in){
  
  r = dim(MCMC_sample[[1]]$curr_param$M)[1]
  c = dim(MCMC_sample[[1]]$curr_param$M)[2]
  nclust = dim(MCMC_sample[[1]]$curr_param$M)[3]
  
  F = array(rep(0,r*c*nclust), c(r, c, nclust))
  D = array(rep(0,c*c*nclust), c(c, c, nclust))
  pi_vec = rep(0.0,nclust)
  
  MCMC_len = length(MCMC_sample)
  
  DIC_vec = rep(0.0,MCMC_len-burn_in+1)
  

  
  for(i in burn_in:MCMC_len){
    
    for(k in 1:nclust){
      F[,,k] = MCMC_sample[[i]]$curr_param$M[,,k]%*%MCMC_sample[[i]]$curr_param$D[,,k]%*%MCMC_sample[[i]]$curr_param$V[,,k]
      D[,,k] = MCMC_sample[[i]]$curr_param$D[,,k]
      
    }
    pi_vec = MCMC_sample[[i]]$pi_vec
    
    DIC_vec[i-burn_in+1] = -2*calculate_log_like_for_data(pi_vec,F,D,data,nclust)
    #print(paste0("DIC_vec[i] = ",DIC_vec[i]))
  }
  
  D_bar = mean(DIC_vec)
  D_var = var(DIC_vec)
  
  DIC = D_bar + D_var/2
  
  return(DIC)
}







calculate_log_like_for_data <-function(pi_vec,F,D,data,nclust){
  data_len = dim(data)[3]
  log_like_for_one_mcmc_iter = 0.0
  
  dyn.load("./C_code/functions_hyper_2by2_R.so")
  hyper0F1_val_arr = rep(0.0,nclust)
  
  for(j in 1:nclust){
    dRet=0.0
    eigenValues = diag(D[,,j])^2/4
    hyper0F1_val =hyper_2by2_R(1.5,eigenValues)
    hyper0F1_val_arr[j] = hyper0F1_val
  }
  
  for(i in 1:data_len){
    log_like_for_one_mcmc_iter = log_like_for_one_mcmc_iter + log(sum(pi_vec * calculate_log_like_for_one_data(F,D,data[,,i],nclust,hyper0F1_val_arr)))
  }
  
  return(log_like_for_one_mcmc_iter)
}





calculate_log_like_for_data_alt <-function(pi_vec, F, D, data, nclust){
  data_len = dim(data)[3]
  log_like_for_one_mcmc_iter = 0.0
  
  
  hyper0F1_val_arr = rep(0.0,nclust)
  
  for(j in 1:nclust){
    dRet=0.0
    eigenValues = diag(D[,,j])^2/4
    hyper0F1_val = hyper_2by2_R(1.5,eigenValues)
    hyper0F1_val_arr[j] = hyper0F1_val
  }
  
  for(i in 1:data_len){
    log_like_for_one_mcmc_iter = log_like_for_one_mcmc_iter + log(sum(pi_vec * calculate_log_like_for_one_data(F,D,data[,,i],nclust,hyper0F1_val_arr)))
  }
  
  return(log_like_for_one_mcmc_iter)
}







calculate_log_like_for_one_data <- function(F,D,data_i,nclust,hyper0F1_val_arr){
  
  log_hyper0F1_val_arr = log(hyper0F1_val_arr)
  #browser()
  val_i = rep(0.0,nclust)
  for(j in 1:nclust){
    dRet = 0.0
    eigenValues = (diag(D[,,j]))^2/4
    
    val_i[j] = sum(diag(t(F[,,j])%*%data_i)) - log_hyper0F1_val_arr[j]
    
  }
  exp_val_i = exp(val_i)
  if(sum(is.nan(exp_val_i)) > 0){
    err_msg("calculate_log_like_for_one_data: density calculation error (NaN) !!")
    stop()
  }
  
  if(sum(is.infinite(exp_val_i)) > 0){
    err_msg("calculate_log_like_for_one_data: density calculation error (Inf) !!")
    stop()
  }
  return(exp_val_i)
}




calc_F_dist <-function(MCMC_sample,L,cid){
  
  mcmc_len = length(MCMC_sample)
  f_dist = rep(0.0,mcmc_len)
  
  for(i in 1:mcmc_len){
    M = MCMC_sample[[i]]$curr_param$M[,,cid]
    D = MCMC_sample[[i]]$curr_param$D[,,cid]
    V = MCMC_sample[[i]]$curr_param$V[,,cid]
    
    f_dist[i] = norm(L$curr_param$F[,,cid]-M%*%D%*%V)
  }
  return(f_dist)
}





calc_F_avg <-function(MCMC_sample,L,cid){
  
  mcmc_len = length(MCMC_sample)
  d = dim(MCMC_sample[[1]]$curr_param$M)
  f_avg = matrix(rep(0.0,d[1]*d[2]),ncol=d[2])
  
  for(i in 1:mcmc_len){
    M = MCMC_sample[[i]]$curr_param$M[,,cid]
    D = MCMC_sample[[i]]$curr_param$D[,,cid]
    V = MCMC_sample[[i]]$curr_param$V[,,cid]
    
    F = M%*%D%*%V
    f_avg = f_avg+F
  }
  return(f_avg/mcmc_len)
}



