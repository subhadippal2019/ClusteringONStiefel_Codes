calc_mis_class_cnt =  function(true_arr,est_arr){
  #### need to adjust for label switching
  u = unique(true_arr)
  
  library(combinat)
  list_per = permn(u)
  k = length(list_per)
  min_mis_cnt = length(true_arr)
  for(i in 1:k){
    mis_cnt = sum(list_per[[i]][true_arr] != est_arr)
    if(mis_cnt < min_mis_cnt)
      min_mis_cnt = mis_cnt
  }
  return(min_mis_cnt)
}

simulated_data_result_summary = function(true_n_cluster,guess_n_cluster,data_id = 1){
  #### loading L
  fName_L = sprintf("~/Desktop/DTI_project/simulated_data/ML_dataset_%d_%d.RData",true_n_cluster,data_id)
  load(fName_L)
  #### loading output file
  fName_output = sprintf("~/Desktop/DTI_project/simulated_data/results_nclust_%d/outputs/output_nclust_%d_data_id_%d_true_nClust_%d.RData",
            true_n_cluster,guess_n_cluster,data_id,true_n_cluster)
  load(fName_output)
  #print(fName)
  
  # total number of data 
  n_data = length(L$clust)
  u = unique(L$clust)
  mis_class_cnt = 0
  burn_in = 1000
  mcmc_len = length(data_with_init_with_MCMC_samples$MCMC_sample)
  for (k in (burn_in+1):mcmc_len){
    mis_class_cnt = mis_class_cnt + calc_mis_class_cnt(L$clust,data_with_init_with_MCMC_samples$MCMC_sample[[k]]$curr_param$id_arr)
  }
  used_mcmc_len = mcmc_len-burn_in
  mis_class_rate = mis_class_cnt / (used_mcmc_len*n_data)
    
  return(mis_class_rate)  
}


mis_class_3 = rep(0,50)
mis_class_4 = rep(0,50)

true_n_cluster = c(3,4)

for(i in true_n_cluster){
  for(data_id in 1:50){
    if(i == 3){
      load("DIC_table_for_nclust_3.RData")
      arr = apply(DIC_table,1,mean)
      j = 1+which.min(arr)
      mis_class_3[data_id] = simulated_data_result_summary(i,j,data_id)
    }
    if(i == 4){
      load("DIC_table_for_nclust_4.RData")
      arr = apply(DIC_table,1,mean)
      j = 1+which.min(arr)
      mis_class_4[data_id] = simulated_data_result_summary(i,j,data_id)
    }
  }### data_id
}## i
  
cat("avg. mis-classification true_n_clust(3) is ",mean(mis_class_3)," so on an avg.",ceiling(mean(mis_class_3)*400)," errors out of 400","\n")
cat("avg. mis-classification true_n_clust(4) is ",mean(mis_class_4)," so on an avg.",ceiling(mean(mis_class_4)*500)," errors out of 500","\n")

load("DIC_table_for_nclust_3.RData")
apply(DIC_table,2,function(x) which.min(x)) 
sum(apply(DIC_table,2,function(x) which.min(x)) != 2)

load("DIC_table_for_nclust_4.RData")
apply(DIC_table,2,function(x) which.min(x)) 
sum(apply(DIC_table,2,function(x) which.min(x)) != 3)

