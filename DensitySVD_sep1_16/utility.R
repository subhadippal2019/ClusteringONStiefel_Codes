###########
### utility functions

print_debug <-function(str,debug_flag){
  if(debug_flag == 1){
    print(str)
  }
}

err_msg <-function(str){
  print(paste0('***********ERROR*************'))
  print(str)
  print(paste0('***********ERROR*************'))
}  


load_src_libs <-function()
{
  library(rstiefel)
  library(gtools)
  library(expm)
  library(MASS)
  source("cluster_param_update.R")
  source("rdensity_d1_d2.R")
  source("NR_method_d1.R")
  source("analyse_MCMC_sample.R")
  source("bivariate_NR_method_d1_d2.R")
  source("stiefel_SVD.R")
  source("utility.R")
  source("preprocess_data.R")
  source("generate_simulated_data_ML.R")
  source("finiteMixtureML.R")
  
  dyn.load("./C_code/functions_hyper_2by2_R.so")
}


### large \phi approximation from Mardia/Khatri 1977 equation 2.19
approx_hyper <-function(eigenvalues){
  
  D = 2*sqrt(eigenvalues)
  
  phi1 = D[1]
  phi2 = D[2]
  
  const_term = (-0.5)*log(2) - log(pi) + log(gamma(1.5))
  out = const_term + (phi1+phi2) - 0.5 * log(phi1+phi2) - 0.5 * (log(phi1) + log(phi2))
  
  return(exp(out))
  
}

approx_hyper_log <-function(eigenvalues){
  
  D = 2*sqrt(eigenvalues)
  
  phi1 = D[1]
  phi2 = D[2]
  
  const_term = (-0.5)*log(2) - log(pi) + log(gamma(1.5))
  out = const_term + (phi1+phi2) - 0.5 * log(phi1+phi2) - 0.5 * (log(phi1) + log(phi2))
  
  return(out)
  
}

test_approx <- function(){
  V1 = NULL
  V2 = NULL
  for(i in 1:5){
    D = c(70+i,60+i)
    eigenvalues = D^2/4
    hyper0F1_val = .C("hyper_2by2_R",cc=1.5,eigenvalues,dRet=0.0)[[3]]
    V1[i] = log(hyper0F1_val)
    V2[i] = log(approx_hyper(eigenvalues))
  }
  print(V1-V2)
  plot(V1,V2)
}

calculate_soft_membership = function(data_with_init_with_MCMC_samples,mcmc_iter_id){

  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample 
  
  data = data_with_init_with_MCMC_samples$data
  curr_param = MCMC_sample[[mcmc_iter_id]]$curr_param
  pi_vec = MCMC_sample[[mcmc_iter_id]]$pi_vec
  N = length(MCMC_sample[[mcmc_iter_id]]$curr_param$id_arr)
  nc = length(pi_vec)
  d_tmp_val = matrix(rep(0,N*nc),ncol=nc)
  
  soft_cluster_for_data = matrix(rep(0,N*nc),nrow=N)
  for (cluster_id in 1:nc){
    d_tmp_val[,cluster_id] = apply(data,3,function(x) dMatrixLangevin(x,curr_param$M[,,cluster_id],curr_param$D[,,cluster_id],curr_param$V[,,cluster_id]))
  }    
  soft_cluster_for_data = d_tmp_val*t(matrix(rep(pi_vec,N),ncol=N))
  soft_cluster_for_data = soft_cluster_for_data/rowSums(soft_cluster_for_data)
  return(soft_cluster_for_data)
}


#sim4=read.table("~/Dropbox/ClusteringDTIonstiefel/DensitySVD_sep1_16/sim_true_4_dic_4_5_6_7_8.txt")
#DIC5_sim4 = sim4[,4]
#DIC5_sim3_mat_50_4 = matrix(DIC5_sim3,nrow=50)
#DIC5_sim4_mat_50_5_new=DIC5_sim4_mat_50_5
#for(i in 1:5) {DIC5_sim4_mat_50_5_new[,i] =  DIC5_sim4_mat_50_5[,i] + (i+1)*6*10*log(500)}
#table(apply(DIC5_sim4_mat_50_5_new,1,function(x) which.min(x)))