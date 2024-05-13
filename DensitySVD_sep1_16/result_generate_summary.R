result_generate_summary <-function(true_nclust = 3){
  
  source("utility.R")
  load_src_libs()
  
  nid = 2
  if(true_nclust == 3){
    n_clust = c(2,3,4,5)
  }
  if(true_nclust == 4){
    n_clust = c(2,3,4,5,6)
  }
  DIC_table = matrix(rep(0,nid*length(n_clust)),ncol=nid)
  
  burn_in = 500
  #output_nclust_3_data_id_6_true_nClust_3.RData 
  for (i in 1:nid){
    for(j in n_clust){
      fName = sprintf("~/Desktop/results_nclust_%d/outputs/output_nclust_%d_data_id_%d_true_nClust_%d.RData",true_nclust,j,i,true_nclust)
      load(fName)
      MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
      data = data_with_init_with_MCMC_samples$data
      
      print(paste0("i = ",i," j = ",j))
      k = j-1
      DIC_table[k,i] = calculate_DIC(MCMC_sample,data,burn_in)
    }       
  }
  save(DIC_table, file="DIC_table_for_nclust_3.RData")
}

result_generate_summary_real <-function(id_str = 12345){
  
  source("utility.R")
  load_src_libs()
  
  n_clust = c(6,7,8)
  DIC_table = rep(0,length(n_clust))
  
  burn_in = 995
  #load("output_nclust_6_real_data_id_12346.RData")
  for(j in 1:length(n_clust)){
    fName = sprintf("~/Desktop/results_real_6_7_8_1000/outputs/output_nclust_%d_real_data_id_%s.RData",n_clust[j],id_str)
    print(fName)
    load(fName)
    MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
    data = data_with_init_with_MCMC_samples$data
    
    print(paste0("j = ",j))

    DIC_table[j] = calculate_DIC(MCMC_sample,data,burn_in)
  }       
  
  file_name = sprintf("DIC_table_for_%s.RData",id_str)
  save(DIC_table, file=file_name)
}


F_compare_summary <- function(true_nclust = 3){
  Nid =  50
  F_est = rep(0.0,Nid)
  F_dist = NULL
  burn_in = 500
  for (i in 1:Nid){
    fName2 = sprintf("~/Desktop/simulated_data/ML_dataset_%d_%d.RData",true_nclust,i)
    load(fName2)
    F_true = L$curr_param$F
    
    fName1 = sprintf("~/Desktop/results_nclust_%d/outputs/output_nclust_%d_data_id_%d_true_nClust_%d.RData",true_nclust,true_nclust,i,true_nclust)
    load(fName1)
    MCMC_sample_0 = data_with_init_with_MCMC_samples$MCMC_sample
    MCMC_sample = MCMC_sample_0[burn_in:length(MCMC_sample_0)]
    l1 = lapply(MCMC_sample,function(x){find_permutation(F_true,x$curr_param$F,true_nclust)})
    M = t(matrix((unlist(l1)),nrow=true_nclust))
    k = 1:length(MCMC_sample)
    M_c = cbind(M,k)
    l2 = lapply(apply(M_c, 1, function(x) list(x)), function(y) unlist(y))
    
    F_list = lapply(l2,function(x){MCMC_sample[[x[(true_nclust+1)]]]$curr_param$F[,,x[1:true_nclust]]})
    
    #F_list = (lapply(MCMC_sample,function(x){x$curr_param$F}))
    F_est = Reduce("+",F_list)/length(F_list)
    
    #p_vec = find_permutation(F_true,F_est)
    #print(paste0("i = ",i))
    #print(p_vec)
    F_dist[i] = Reduce("+",(F_true - F_est)^2)
  }
  return(F_dist)
}
  
  find_permutation <- function(F_true,F_est,true_nclust){
    F1 = F_true
    F2 = F_est
    perm_vector = 1:true_nclust
    p_list = permn(perm_vector)
    min_dist = 99999
    min_dist_idx = 1
    for (i in 1:length(p_list)){
      F_tmp = F2[,,p_list[[i]]]
      d1=Reduce("+",(F1 - F_tmp)^2)
      if(d1 < min_dist){
        min_dist = d1
        min_dist_idx = i 
      }
    }
    return(p_list[[min_dist_idx]])
  }


