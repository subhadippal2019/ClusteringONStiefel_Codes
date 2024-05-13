##### MLE check

generate_data_ML<-function(N=10,data_dir="./testMLE"){
  
  source("bivariate_NR_method_d1_d2.R")
  source("stiefel_SVD.R")
  ###################
  ### Description:
  
  ### Simulated data generation for fixed number of clusters 
  ### Visualize in some random direction on 2D
  ### Clustering heatmap visualization
  
  ###################
  
  ##################
  library(rstiefel)
  ##################
  
  
  nClust = 1
  
  M_hyper_arr = array(c(1,0,0,0,0,1), c(3, 2, nClust))
  
  ### mean matrix generation nClust matrices
  M = array(rep(0,3*2*nClust), c(3, 2, nClust))
  for(i in 1:nClust){
    #M[,,i] = rmf.matrix(M_hyper_arr[,,i]) 
    M[,,i] = M_hyper_arr[,,i]
  }
  
  D = array(rep(0,2*2*nClust), c(2, 2, nClust))
  for(i in 1:nClust){
    s = 0.5*sort(rgamma(c(1,1),c(160,100),c(20,20)),decreasing=TRUE)
    D[,,i] = diag(s+i)
  }
  
  V = array(rep(0,2*2*nClust), c(2, 2, nClust))
  for(i in 1:nClust){
    V[,,i] = matrix(c(1,0,0,1))
  }
  ###### main parameter
  F = array(rep(0,3*2*nClust), c(3, 2, nClust))
  
  ##########################
  ### to keep the sign thing consistent
  
  for(i in 1:nClust){
    M_V_list = change_sign_M_V(M[,,i],V[,,i])
    M[,,i] = M_V_list$M
    V[,,i] = M_V_list$V
    F[,,i] = M[,,i]%*%D[,,i]%*%V[,,i]
  }
  
  
  ##########################
  data = array(rep(0,3*2*N), c(3, 2, N)); 
  ## 0.2 0.4 0.4
  for(i in 1:N){
    F1 = M[,,1]%*%D[,,1]%*%V[,,1]
    
    data[,,i] = rmf.matrix(F1) 
  }
  
  #####
  
  data_param = list(data=data,M=M,D=D,V=V)
  fileName = sprintf("%s/data_test_%d.rda",data_dir,N)
  save(data_param,file=fileName)
  
}
  
  
initial_parameter_est <-function(data){
  
  library(expm)
  return_all_cluster_initial_est = NULL
  
  # sum_X = matrix(rep(0.0,3*2),c(3,2))
  X_bar = matrix(rep(0.0,3*2),c(3,2))
  #   L1 = length(idx)
  #   for(i in 1:L1){
  #     sum_X = sum_X + L$data[,,idx[i]]
  #   }
  #   X_bar = sum_X/L1
  #   
  X_bar = apply(data,c(1,2),mean)
  
  #   #%%%% following Mardia's notation
  #   R_bar = sqrtm(t(X_bar)%*%X_bar)
  #   polar_M = X_bar%*%solve(R_bar)
  #   eig = eigen(R_bar)  ##### R_bar = t(U)*D*U  
  #   U=t(eig$vectors)
  #   g = eig$values ## as a vector
  #   M = polar_M%*%t(U)
  #   
  ### following Chikuse book page 111 (for M and V as well)
  res = stiefel_SVD(X_bar) ### X_bar = M*D*V
  
  M = res$u
  return_all_cluster_initial_est$M = M
  
  V = res$v
  return_all_cluster_initial_est$V = V
  
  #### phi needs to be solved
  ### for small phi, phi = p*g_i ************************** TO BE DONE 
  p = dim(X_bar)[2]
  
  ##if(phi are small){
  #phi = p*g 
  #}
  phi = bivariate_NR_method_d1_d2(p,res$d[1],res$d[2],max_iter=50)
  
  ##if(phi are small){
  #phi = 2(1-g)
  #}
  cat("phi:",phi,"\n")
  
  cat("res$d:",res$d[1],res$d[2],"\n")
  
  D = diag(phi)
  return_all_cluster_initial_est$D = D
  
  
  return(return_all_cluster_initial_est) ### ***** NOT WORKING WELL
}

#### main code
source("stiefel_SVD.R")
source("bivariate_NR_method_d1_d2.R")
dyn.load("./C_code/functions_hyper_2by2_R.so")
KK = 100

dist_M = rep(0,KK)
dist_D = rep(0,KK)
dist_V = rep(0,KK)

fileName = sprintf("./testMLE/data_test_10000.rda")
load(fileName)

for (i in 1:KK){
  
  res = initial_parameter_est(data_param$data[,,(1:(i*100))])
  
  est_M = res$M
  est_D = res$D
  est_V = res$V
  
  orig_M = data_param$M[,,1]
  orig_D = data_param$D[,,1]
  orig_V = data_param$V[,,1]
  
  
  p = dim(orig_M)[2]  
  
  dist_M[i] = p - sum(diag(t(orig_M)%*%est_M))
  
  dist_V[i] = p - sum(diag(t(orig_V)%*%est_V))
  
  tmp = diag(orig_D) - diag(est_D)
  
  
  dist_D[i] = sqrt(t(tmp)%*%tmp)
  
  
  
  cat("done ",i,"\n")
}




