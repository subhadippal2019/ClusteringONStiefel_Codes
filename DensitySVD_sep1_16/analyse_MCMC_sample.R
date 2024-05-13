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
  
  #l1 = lapply(MCMC_sample,function(x){(F_true,x$curr_param$F,true_nclust)})
  
  MCMC_sample_1 = MCMC_sample[burn_in:MCMC_len]
  DIC_vec = lapply(MCMC_sample_1,function(x){-2*calculate_log_like_for_data(x$pi_vec,x$curr_param$F,x$curr_param$D,data,nclust)} )
  #DIC_vec = lapply(MCMC_sample_1,function(x){-2*test1(x$curr_param$F,data)} )
  DIC_vec = unlist(DIC_vec)
 
  #for(i in burn_in:MCMC_len){
    
    #for(k in 1:nclust){
      ##F[,,k] = MCMC_sample[[i]]$curr_param$M[,,k]%*%MCMC_sample[[i]]$curr_param$D[,,k]%*%MCMC_sample[[i]]$curr_param$V[,,k]
      #F[,,k] = MCMC_sample[[i]]$curr_param$F[,,k]
      #D[,,k] = MCMC_sample[[i]]$curr_param$D[,,k]
      
    #}
    #F = MCMC_sample[[i]]$curr_param$F
    #D = MCMC_sample[[i]]$curr_param$D
    
    #pi_vec = MCMC_sample[[i]]$pi_vec
    
    #DIC_vec[i-burn_in+1] = -2*calculate_log_like_for_data(pi_vec,F,D,data,nclust)
    #print(paste0("DIC_vec[i] = ",DIC_vec[i]))
  #}
  
  D_bar = mean(DIC_vec)
  D_var = var(DIC_vec)
  
  DIC = D_bar + D_var/2
  
  return(DIC)
}

calculate_DIC_with_penalty <-function(MCMC_sample,data,burn_in){
  
  r = dim(MCMC_sample[[1]]$curr_param$M)[1]
  c = dim(MCMC_sample[[1]]$curr_param$M)[2]
  nclust = dim(MCMC_sample[[1]]$curr_param$M)[3]
  
  F = array(rep(0,r*c*nclust), c(r, c, nclust))
  D = array(rep(0,c*c*nclust), c(c, c, nclust))
  pi_vec = rep(0.0,nclust)
  
  MCMC_len = length(MCMC_sample)
  
  DIC_vec = rep(0.0,MCMC_len-burn_in+1)
  
  #l1 = lapply(MCMC_sample,function(x){(F_true,x$curr_param$F,true_nclust)})
  
  MCMC_sample_1 = MCMC_sample[burn_in:MCMC_len]
  DIC_vec = lapply(MCMC_sample_1,function(x){-2*calculate_log_like_for_data(x$pi_vec,x$curr_param$F,x$curr_param$D,data,nclust)} )
  #DIC_vec = lapply(MCMC_sample_1,function(x){-2*test1(x$curr_param$F,data)} )
  DIC_vec = unlist(DIC_vec)
  
  #for(i in burn_in:MCMC_len){
  
  #for(k in 1:nclust){
  ##F[,,k] = MCMC_sample[[i]]$curr_param$M[,,k]%*%MCMC_sample[[i]]$curr_param$D[,,k]%*%MCMC_sample[[i]]$curr_param$V[,,k]
  #F[,,k] = MCMC_sample[[i]]$curr_param$F[,,k]
  #D[,,k] = MCMC_sample[[i]]$curr_param$D[,,k]
  
  #}
  #F = MCMC_sample[[i]]$curr_param$F
  #D = MCMC_sample[[i]]$curr_param$D
  
  #pi_vec = MCMC_sample[[i]]$pi_vec
  
  #DIC_vec[i-burn_in+1] = -2*calculate_log_like_for_data(pi_vec,F,D,data,nclust)
  #print(paste0("DIC_vec[i] = ",DIC_vec[i]))
  #}
  
  D_bar = mean(DIC_vec)
  D_var = var(DIC_vec)
  
  penalty = 0
  DIC = D_bar + D_var/2 + penalty
  
  return(DIC)
}


calculate_log_like_for_data <-function(pi_vec,F,D,data,nclust){
  data_len = dim(data)[3]
  log_like_for_one_mcmc_iter = 0.0
  
  #dyn.load("./C_code/functions_hyper_2by2_R.so")
  hyper0F1_val_arr = rep(0.0,nclust)
  
  #for(j in 1:nclust){
  #  dRet=0.0
  #  eigenValues = diag(D[,,j])^2/4
  #  hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,dRet)[[3]]
  #  hyper0F1_val_arr[j] = hyper0F1_val
  #}
  dRet = 0.0
  hyper0F1_val_arr = apply(D,3,function(x){.C("hyper_2by2_R",a=1.5,diag(x)^2/4,dRet)[[3]]})
  
  vec_data_log_like = apply(data,3,function(x) {log(sum(pi_vec*calculate_like_for_one_data(F,D,x,nclust,hyper0F1_val_arr)))} )
  log_like_for_one_mcmc_iter = sum(vec_data_log_like)
  #for(i in 1:data_len){
  #  log_like_for_one_mcmc_iter = log_like_for_one_mcmc_iter + log(sum(pi_vec * calculate_like_for_one_data(F,D,data[,,i],nclust,hyper0F1_val_arr) ))
  #}
  
  return(log_like_for_one_mcmc_iter)
}


calculate_like_for_one_data <- function(F,D,data_i,nclust,hyper0F1_val_arr){
  
  log_hyper0F1_val_arr = log(hyper0F1_val_arr)
  #browser()
  val_i = rep(0.0,nclust)
  #for(j in 1:nclust){
    ##dRet = 0.0
    ##eigenValues = (diag(D[,,j]))^2/4
    
    #val_i[j] = sum(diag(t(F[,,j])%*%data_i)) - log_hyper0F1_val_arr[j]
    
  #}
  val_i = apply(F,3,function(x){sum(diag(t(x)%*%data_i))} ) - log_hyper0F1_val_arr
  
  exp_val_i = exp(val_i)
  if(sum(is.nan(exp_val_i)) > 0){
    err_msg("calculate_like_for_one_data: density calculation error (NaN) !!")
    stop()
  }
  
  if(sum(is.infinite(exp_val_i)) > 0){
    err_msg("calculate_like_for_one_data: density calculation error (Inf) !!")
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

calculate_z_hat <- function(MCMC_sample){
  mcmc_len = length(MCMC_sample)
  N = length(MCMC_sample[[1]]$curr_param$id_arr)
  z_hat_mat = matrix(rep(0,N*mcmc_len),ncol=mcmc_len)
  for(i in 1:mcmc_len){
    z_hat_mat[,i] = MCMC_sample[[i]]$curr_param$id_arr
  }
  all_data_z = rep(0,N)
  all_data_z = apply(z_hat_mat,1,function(x) find_z_mode(x))
  
  return(all_data_z)
}

find_z_mode <-function(Z){
  return(as.integer(names(which.max(table(Z)))))
}

calculate_theta_hat_all_cluster <- function(MCMC_sample){
  r = dim(MCMC_sample[[1]]$curr_param$M)[1]
  c = dim(MCMC_sample[[1]]$curr_param$M)[2]
  nclust = dim(MCMC_sample[[1]]$curr_param$F)[3]
  F_hat_arr = array(rep(0,r*c*nclust), c(r, c, nclust))
  
  for(cid in 1:nclust){ #### TBC 
    F_hat_arr[,,cid] = calculate_theta_hat(MCMC_sample,cid)
  }
  return(F_hat_arr)  
}

calculate_theta_bar_all_cluster <- function(MCMC_sample){
  r = dim(MCMC_sample[[1]]$curr_param$M)[1]
  c = dim(MCMC_sample[[1]]$curr_param$M)[2]
  
  
  nclust = dim(MCMC_sample[[1]]$curr_param$F)[3]
  F_bar_arr = array(rep(0,r*c*nclust), c(r, c, nclust))
  
  for(cid in 1:nclust){
    F_bar_arr[,,cid] = calculate_theta_bar(MCMC_sample,cid)
  }
  return(F_bar_arr)  
}

calculate_theta_hat <- function(MCMC_sample,cid){
  mcmc_len = length(MCMC_sample)
  r = dim(MCMC_sample[[1]]$curr_param$M)[1]
  c = dim(MCMC_sample[[1]]$curr_param$M)[2]
  
  F_arr = array(rep(0,r*c*mcmc_len), c(r, c, mcmc_len))
  for(i in 1:mcmc_len){
    F_arr[,,i] = MCMC_sample[[i]]$curr_param$F[,,cid]
  }
  F_hat = apply(F_arr,c(1,2),function(x) compute_mode(x))
  return(F_hat)
}

calculate_theta_bar <- function(MCMC_sample,cid){
  mcmc_len = length(MCMC_sample)
  r = dim(MCMC_sample[[1]]$curr_param$M)[1]
  c = dim(MCMC_sample[[1]]$curr_param$M)[2]
  
  F_arr = array(rep(0,r*c*mcmc_len), c(r, c, mcmc_len))
  for(i in 1:mcmc_len){
    F_arr[,,i] = MCMC_sample[[i]]$curr_param$F[,,cid]
  }
  F_bar = apply(F_arr,c(1,2),mean)
  return(F_bar)
}

compute_mode <-function(X){
  G = density(X)
  return(G$x[which.max(G$y)])
}

MCMC_F_hat_trace <-function(nclust,dataset_id,clust_id,element_id){
  burn_in = 1
  fileName1 = sprintf("~/Desktop/DTI_project/simulated_data/ML_dataset_%d_%d.RData",nclust,dataset_id)
  load(fileName1)
  F_t = L$curr_param$F[,,clust_id]
  F_t_vec = matrix(F_t,nrow=6)
  
  
  fileName2 = sprintf("~/Desktop/DTI_project/simulated_data/results_nclust_%d/outputs/output_nclust_%d_data_id_%d_true_nClust_%d.RData",nclust,nclust,dataset_id,nclust)
  load(fileName2)
  mcmc_len = length(data_with_init_with_MCMC_samples$MCMC_sample)
  mcmc_sample = data_with_init_with_MCMC_samples$MCMC_sample[burn_in:mcmc_len]
  mcmc_len_after_burn_in = length(mcmc_sample)
  arr_mcmc_elem = rep(0,mcmc_len_after_burn_in)
  for(i in 1:mcmc_len_after_burn_in){
    F_vec = matrix(data_with_init_with_MCMC_samples$MCMC_sample[[burn_in+i-1]]$curr_param$F[,,clust_id],nrow=6)
    arr_mcmc_elem[i] = F_vec[element_id]
  }
  arr_mcmc_elem_cumsum = cumsum(arr_mcmc_elem)/(1:length(arr_mcmc_elem))
  fig_name = sprintf("true_%d_data_%d_clust_id_%d_elem_id_%d_acf.pdf",nclust,dataset_id,clust_id,element_id)
  pdf(fig_name)
  arr_mcmc_elem_acf = acf(arr_mcmc_elem,main="",xlab="lag",ylab="Autocorrelation")
  dev.off()
  
  fig_name = sprintf("true_%d_data_%d_clust_id_%d_elem_id_%d_cumsum.pdf",nclust,dataset_id,clust_id,element_id)
  pdf(fig_name)
  plot(arr_mcmc_elem_cumsum,xlab="MCMC iteration",ylab="Cumulative average",type='l')
  abline(F_t_vec[element_id],0,col='blue')
  dev.off()
}


calc_MSE_F__F_hat <- function(){
  library(combinat)
  data_set_len = 50
  MSE_perm_matrix = NULL
  MSE_matrix = matrix(rep(0,2*data_set_len),ncol=data_set_len)
  Perm_matrix = matrix(rep(0,2*data_set_len),ncol=data_set_len)
  burn_in = 500
  for(t_clust in 3:4){
    for(i in 1:data_set_len){
      fileName1 = sprintf("~/Desktop/DTI_project/simulated_data/ML_dataset_%d_%d.RData",t_clust,i)
      load(fileName1)
      F_t = L$curr_param$F
      fileName2 = sprintf("~/Desktop/DTI_project/simulated_data/results_nclust_%d/outputs/output_nclust_%d_data_id_%d_true_nClust_%d.RData",t_clust,t_clust,i,t_clust)
      load(fileName2)
      F_h = array(rep(0,3*2*t_clust), c(3, 2, t_clust))
      for(j in 1:t_clust){
        mcmc_len = length(data_with_init_with_MCMC_samples$MCMC_sample)
        mcmc_sample = data_with_init_with_MCMC_samples$MCMC_sample[burn_in:mcmc_len]
        F_h[,,j] = calculate_theta_bar(mcmc_sample,j)
      }
      
      ret_MSE = MSE_F_t_F_h(F_t,F_h)
      MSE_matrix[t_clust-2,i] = ret_MSE$val
      Perm_matrix[t_clust-2,i] = ret_MSE$perm
    }
  }
  MSE_perm_matrix$MSE_matrix = MSE_matrix
  MSE_perm_matrix$Perm_matrix = Perm_matrix
  
  return(MSE_perm_matrix)
}


MSE_F_t_F_h <-function(F_t,F_h){
  nClust = dim(F_t)[3]
  list_per = permn(nClust)
  K = length(list_per)
  MSE = rep(0,K)
  for(k in 1:K){
    
    F_h_perm = F_h[,,list_per[[k]]]  
    
    S = 0.0
    for(i in 1:nClust){
      S = S+norm(F_t[,,i] - F_h_perm[,,i],type="F")
    }
    S = S/(nClust)
    MSE[k] = S
  }
  
  
  ret_MSE = NULL
  h = which.min(MSE)
  ret_MSE$val = min(MSE)
  ret_MSE$perm = h
  return(ret_MSE)
}

###### alternative DIC calculation
################################# 

calculate_DIC_4 <- function(data_with_init_with_MCMC_samples,burn_in=500){
  data = data_with_init_with_MCMC_samples$data
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  #burn_in = 500
  mcmc_len = length(MCMC_sample)
  MCMC_sample_1 = MCMC_sample[burn_in:mcmc_len]
  
  F_bar_arr = calculate_theta_bar_all_cluster(MCMC_sample_1)
  nclust = dim(F_bar_arr)[3]
  D_bar_arr = array(apply(F_bar_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  tmp_1 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,x$curr_param$F,x$curr_param$D,data,nclust)} )
  
  
  tmp_2 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,F_bar_arr,D_bar_arr,data,nclust)} )
  DIC_val_4 = NULL
  DIC_val_4$val = -4*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  DIC_val_4$pd = -2*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  return(DIC_val_4)
}

calculate_DIC_5 <- function(data_with_init_with_MCMC_samples,burn_in=500){

  data = data_with_init_with_MCMC_samples$data
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  #burn_in = 500
  mcmc_len = length(MCMC_sample)
  MCMC_sample_1 = MCMC_sample[burn_in:mcmc_len]
  
  F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample_1)
  nclust = dim(F_hat_arr)[3]
  D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  tmp_1 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,x$curr_param$F,x$curr_param$D,data,nclust)} )

  Z_hat_arr = calculate_z_hat(MCMC_sample)
  pi_hat_vec = rep(0,nclust)
  T = table(Z_hat_arr)/length(Z_hat_arr)
  pi_hat_vec[as.integer(names(T))]=as.numeric(T)
  tmp_2 = calculate_log_like_for_data_with_z(pi_hat_vec,Z_hat_arr,F_hat_arr,D_hat_arr,data,nclust)
  
  DIC_val_5 = NULL
  DIC_val_5$val = -4*mean(unlist(tmp_1)) + 2*tmp_2
  DIC_val_5$pd = -2*mean(unlist(tmp_1)) + 2*tmp_2
  return(DIC_val_5)
}

calculate_DIC_5_with_penalty <- function(data_with_init_with_MCMC_samples,burn_in=500){
  
  data = data_with_init_with_MCMC_samples$data
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  #burn_in = 500
  mcmc_len = length(MCMC_sample)
  MCMC_sample_1 = MCMC_sample[burn_in:mcmc_len]
  
  F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample_1)
  nclust = dim(F_hat_arr)[3]
  D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  tmp_1 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,x$curr_param$F,x$curr_param$D,data,nclust)} )
  
  Z_hat_arr = calculate_z_hat(MCMC_sample)
  pi_hat_vec = rep(0,nclust)
  T = table(Z_hat_arr)/length(Z_hat_arr)
  pi_hat_vec[as.integer(names(T))]=as.numeric(T)
  tmp_2 = calculate_log_like_for_data_with_z(pi_hat_vec,Z_hat_arr,F_hat_arr,D_hat_arr,data,nclust)
  
  DIC_val_5 = NULL
  DIC_val_5$val = -4*mean(unlist(tmp_1)) + 2*tmp_2
  DIC_val_5$pd = -2*mean(unlist(tmp_1)) + 2*tmp_2
  
  penalty = length(MCMC_sample_1[[1]]$pi_vec) *( dim(MCMC_sample_1[[1]]$curr_param$F)[1]*dim(MCMC_sample_1[[1]]$curr_param$F)[2] )
  
  DIC_val_5_with_penalty = NULL
  DIC_val_5_with_penalty$val = DIC_val_5$val + penalty
  DIC_val_5_with_penalty$pd = DIC_val_5$pd
  return(DIC_val_5_with_penalty)
}


calculate_DIC_6 <- function(data_with_init_with_MCMC_samples,burn_in=500){
  data = data_with_init_with_MCMC_samples$data
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  #burn_in = 500
  mcmc_len = length(MCMC_sample)
  MCMC_sample_1 = MCMC_sample[burn_in:mcmc_len]
  
  F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample_1)
  nclust = dim(F_hat_arr)[3]
  D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  tmp_1 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,x$curr_param$F,x$curr_param$D,data,nclust)} )
  
  tmp_2 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,F_hat_arr,D_hat_arr,data,nclust)} )
  
  DIC_val_6 = NULL
  DIC_val_6$val = -4*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  DIC_val_6$pd = -2*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  
  return(DIC_val_6)
}

calculate_DIC_7 <- function(data_with_init_with_MCMC_samples,burn_in=500){
  data = data_with_init_with_MCMC_samples$data
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  #burn_in = 500
  mcmc_len = length(MCMC_sample)
  MCMC_sample_1 = MCMC_sample[burn_in:mcmc_len]
  
  F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample_1)
  nclust = dim(F_hat_arr)[3]
  D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  tmp_1 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,x$curr_param$F,x$curr_param$D,data,nclust)} )
  
  Z_hat_arr = calculate_z_hat(MCMC_sample)
  pi_hat_vec_all_1 = rep(1,nclust)
  
  tmp_2 = calculate_log_like_for_data_with_z(pi_hat_vec_all_1,Z_hat_arr,F_hat_arr,D_hat_arr,data,nclust)
  
  DIC_val_7 = NULL
  DIC_val_7$val = -4*mean(unlist(tmp_1)) + 2*tmp_2
  DIC_val_7$pd = -2*mean(unlist(tmp_1)) + 2*tmp_2
  return(DIC_val_7)
}

calculate_DIC_8 <- function(data_with_init_with_MCMC_samples,burn_in=500){
  data = data_with_init_with_MCMC_samples$data
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  #burn_in = 500
  mcmc_len = length(MCMC_sample)
  MCMC_sample_1 = MCMC_sample[burn_in:mcmc_len]
  
  F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample_1)
  nclust = dim(F_hat_arr)[3]
  D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  tmp_1 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,x$curr_param$F,x$curr_param$D,data,nclust)} )
  
  pi_hat_vec_all_1 = rep(1,nclust)
  tmp_2 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(pi_vec_all_1,x$curr_param$id_arr,F_hat_arr,D_hat_arr,data,nclust)} )
  
  DIC_val_8 = NULL
  DIC_val_8$val = -4*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  DIC_val_8$pd = -2*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  
  return(DIC_val_8)
}



calculate_DIC_4_5_6_7_8 <-function(data_with_init_with_MCMC_samples,burn_in=500){
  data = data_with_init_with_MCMC_samples$data
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  #burn_in = 500
  mcmc_len = length(MCMC_sample)
  MCMC_sample_1 = MCMC_sample[burn_in:mcmc_len]
  
  F_bar_arr = calculate_theta_bar_all_cluster(MCMC_sample_1)
  nclust = dim(F_bar_arr)[3]
  D_bar_arr = array(apply(F_bar_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  tmp_1 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,x$curr_param$F,x$curr_param$D,data,nclust)} )
  
  tmp_2 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,F_bar_arr,D_bar_arr,data,nclust)} )
  DIC_val = NULL
  DIC_val$val4 = -4*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  DIC_val$pd4 = -2*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample_1)
  nclust = dim(F_hat_arr)[3]
  D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  tmp_2 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(x$pi_vec,x$curr_param$id_arr,F_hat_arr,D_hat_arr,data,nclust)} )
  DIC_val$val6 = -4*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  DIC_val$pd6 = -2*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  
  #F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample)
  #D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  Z_hat_arr = calculate_z_hat(MCMC_sample)
  pi_hat_vec = rep(0,nclust)
  T = table(Z_hat_arr)/length(Z_hat_arr)
  pi_hat_vec[as.integer(names(T))]=as.numeric(T)
  tmp_2 = calculate_log_like_for_data_with_z(pi_hat_vec,Z_hat_arr,F_hat_arr,D_hat_arr,data,nclust)
  
  DIC_val$val5 = -4*mean(unlist(tmp_1)) + 2*tmp_2
  DIC_val$pd5 = -2*mean(unlist(tmp_1)) + 2*tmp_2
  
  #F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample)
  #D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  #Z_hat_arr = calculate_z_hat(MCMC_sample)
  pi_hat_vec_all_1 = rep(1,nclust)
  
  tmp_2 = calculate_log_like_for_data_with_z(pi_hat_vec_all_1,Z_hat_arr,F_hat_arr,D_hat_arr,data,nclust)
  
  DIC_val$val7 = -4*mean(unlist(tmp_1)) + 2*tmp_2
  DIC_val$pd7 = -2*mean(unlist(tmp_1)) + 2*tmp_2
  
  
  #F_hat_arr = calculate_theta_hat_all_cluster(MCMC_sample)
  #D_hat_arr = array(apply(F_hat_arr,3,function(x){diag(sqrt(eigen(t(x)%*%x)$values))}),c(2,2,nclust))
  
  tmp_2 = lapply(MCMC_sample_1,function(x){calculate_log_like_for_data_with_z(pi_hat_vec_all_1,x$curr_param$id_arr,F_hat_arr,D_hat_arr,data,nclust)} )
  DIC_val$val8 = -4*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  DIC_val$pd8 = -2*mean(unlist(tmp_1)) + 2*mean(unlist(tmp_2))
  
  return(DIC_val)
}


##############

calculate_log_like_for_data_with_z <-function(pi_vec,id_arr,F,D,data,nclust){
  #browser()
  data_len = dim(data)[3]
  log_like_for_one_mcmc_iter = 0.0
  
  #dyn.load("./C_code/functions_hyper_2by2_R.so")
  hyper0F1_val_arr = rep(0.0,nclust)
  
  #for(j in 1:nclust){
  #  dRet=0.0
  #  eigenValues = diag(D[,,j])^2/4
  #  hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,dRet)[[3]]
  #  hyper0F1_val_arr[j] = hyper0F1_val
  #}
  dRet = 0.0
  hyper0F1_val_arr = apply(D,3,function(x){.C("hyper_2by2_R",a=1.5,diag(x)^2/4,dRet)[[3]]})
  
  data_with_z = vector("list",data_len)
  for(i in 1:data_len){
    data_with_z[[i]]$data = data[,,i]
    data_with_z[[i]]$z = id_arr[i]
  }
  vec_data_log_like = lapply(data_with_z,function(x) {log(calculate_like_for_one_data_with_z(F,D,x$data,x$z,hyper0F1_val_arr,pi_vec))} )
  log_like_for_one_mcmc_iter = sum(unlist(vec_data_log_like))
  #for(i in 1:data_len){
  #  log_like_for_one_mcmc_iter = log_like_for_one_mcmc_iter + log(sum(pi_vec * calculate_like_for_one_data(F,D,data[,,i],nclust,hyper0F1_val_arr) ))
  #}
  
  return(log_like_for_one_mcmc_iter)
}

calculate_like_for_one_data_with_z <- function(F,D,data_i,z_i,hyper0F1_val_arr,pi_vec){
  F_c = F[,,z_i]
  
  log_hyper0F1_val_arr = log(hyper0F1_val_arr)
  #browser()
  #val_i = rep(0.0,nclust)
  #for(j in 1:nclust){
  ##dRet = 0.0
  ##eigenValues = (diag(D[,,j]))^2/4
  
  #val_i[j] = sum(diag(t(F[,,j])%*%data_i)) - log_hyper0F1_val_arr[j]
  
  #}
  val_i = sum(diag(t(F_c)%*%data_i)) - log_hyper0F1_val_arr[z_i]
  
  exp_val_i = exp(val_i)
  if(sum(is.nan(exp_val_i)) > 0){
    err_msg("calculate_like_for_one_data: density calculation error (NaN) !!")
    stop()
  }
  
  if(sum(is.infinite(exp_val_i)) > 0){
    err_msg("calculate_like_for_one_data: density calculation error (Inf) !!")
    stop()
  }
  return(pi_vec[z_i]*exp_val_i)
}