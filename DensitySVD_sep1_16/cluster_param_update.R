
cluster_param_update<-function(data, curr_param, cluster_assign_vec, cluster_id_vec=c(1),hyper){

  ###################
  ### Description:
  
  ### update of all the parameters one by one 
  ### (first D, then M , then V) for all the clusters
  
  ###################
  
  # setwd("~/Dropbox/ClusteringDTIonstiefel/DensitySVD_sep1_16/")
  # L = generate_simulated_data_ML(N=100) 
  # load("ML_dataset.Rdata")
  # new_param = cluster_param_update(L$data,L$curr_param,L$clust,c(2,3,1),hyper)
  # D=c(5.4,0.7);S=c(14,15);hyper=NULL;hyper$alpha=1;hyper$beta=0
  
  ### cluster_id_vec contains all the id for those clusters that needed to be udpated (e.g: c(3,1,2))
  for(cluster_id in cluster_id_vec){
    
    ID_cluster = which(cluster_assign_vec == cluster_id)
    N = length(ID_cluster)
    if(N == 0){
      sum_X = data[,,1]*0
    }else{
      X = data[,,ID_cluster]
      sum_X = apply(X,c(1,2),sum)
#       t1= Sys.time()
#       if(N == 1){
#         sum_X = X
#       }else{
#         sum_X = X[,,1]*0
#         for (i in 1:N){
#           sum_X = sum_X + X[,,i]
#         }
#       }
#       Sys.time()-t1
    }
    ####################################
    ############# inference ############
    ####################################
  
    M = curr_param$M[,,cluster_id]
    D = curr_param$D[,,cluster_id]
    V = curr_param$V[,,cluster_id]
    
    if(N > 5){ #### if cluster has at least one data point; how many points?? 
               #### for practical purpose may be we can take N>5 to avoid the case S being close to N
      ####################################
      
      ### old M,D,V parameter update
      ####################################
      
      #### conditional of D|M,V for Gibbs
      S = t(M)%*%sum_X%*%t(V)
      ### for D[1] = d_1
      #print(D)
      #print(S)
      
      #####%%%%%%%%%#########**********************
      #### here according to cluster size parameter beta is changing
      hyper$beta = N*0.004
      #############################################
      d1 = rdensity_d1_d2(D,S,dimIndex=1,N,nBin=100,sampleSize=1,hyper,cluster_id)
      D[1,1] = d1
      print_debug(paste0("new d1 = ",d1),hyper$debug)
      ### for D[2] = d_2
      d2 = rdensity_d1_d2(D,S,dimIndex=2,N,nBin=100,sampleSize=1,hyper,cluster_id)
      D[2,2] = d2
      print_debug(paste0("new d2 = ",d2),hyper$debug)
      
      #### conditional of M|D,V for Gibbs
      
      M_param = hyper$G + sum_X%*%t(V)%*%(D)
      #print(M_param)
      M = rmf.matrix(M_param) 
      ### TBD: make the max absolute value for columnwise vector positive 
      
      #### this should automatically make V with proper sign
      M = change_sign_M(M)
     
      #### conditional of V|M,D for Gibbs
      V_param = hyper$H + (D)%*%t(M)%*%sum_X
      V = rmf.matrix(V_param) 
      ### TBD: make the max absolute value for columnwise vector positive 
    } ### if N > 0
    
    #### update the current parameter set
    
    curr_param$M[,,cluster_id] = M
    curr_param$D[,,cluster_id] = D
    curr_param$V[,,cluster_id] = V
    
    curr_param$F[,,cluster_id] = M%*%D%*%V
    curr_param$id_arr = cluster_assign_vec
    print(paste0("cluster update done for ",cluster_id))
  }
  
  return(curr_param)
  

}
