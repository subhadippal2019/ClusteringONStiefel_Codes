#### VCG dataset
#### http://www.mcs.st-and.ac.uk/~pej/vcg.html
source("utility.R")
load_src_libs()

### find q and p from Frank system
### 


vcg_data = read.csv("./vcg_data/vcg.csv")

X1 = vcg_data[,3:5]
X2 = vcg_data[,6:8]

X1_M = vcg_data[,12:14]
X2_M = vcg_data[,15:17]

##### Entire data
output_file = "./vcg_data/vcg_output.RData"

N = 98
Y = array(c(0,0,0,0,0,0), c(3, 2, N))
Y_M = Y

for (i in 1:N){
  Y[,1,i] = t(X1[i,])
  Y[,2,i] = t(X2[i,])
  
  Y_M[,1,i] = t(X1_M[i,])
  Y_M[,2,i] = t(X2_M[i,])
}

data = Y
data_M = Y_M
#data_with_init_with_MCMC_samples = finiteMixtureML(data,nc=1,true_nc=1,max_iter=5000,run_id=1,output_file,vague_prior=3)

#########

### for each of the four age groups
grp_1_id = which(vcg_data$AgeSex==1)
grp_2_id = which(vcg_data$AgeSex==2)
grp_3_id = which(vcg_data$AgeSex==3)
grp_4_id = which(vcg_data$AgeSex==4)


max_iter = 5000

##### group 1 data

X1_1 = vcg_data[grp_1_id,3:5]
X2_1 = vcg_data[grp_1_id,6:8]

X1_1_M = vcg_data[grp_1_id,12:14]
X2_1_M = vcg_data[grp_1_id,15:17]

N = length(grp_1_id)
Y = array(c(0,0,0,0,0,0), c(3, 2, N))
Y_M = Y

for (i in 1:N){
  Y[,1,i] = t(X1_1[i,])
  Y[,2,i] = t(X2_1[i,])
  Y_M[,1,i] = t(X1_1_M[i,])
  Y_M[,2,i] = t(X2_1_M[i,])
}

data = Y
data_M = Y_M

output_file = "./vcg_data/vcg_output_group1.RData"
data_with_init_with_MCMC_samples = finiteMixtureML(data,nc=1,true_nc=1,max_iter,run_id=1,output_file,vague_prior=3)

##### group 2 data
X1_2 = vcg_data[grp_2_id,3:5]
X2_2 = vcg_data[grp_2_id,6:8]

X1_2_M = vcg_data[grp_2_id,12:14]
X2_2_M = vcg_data[grp_2_id,15:17]


N = length(grp_2_id)
Y = array(c(0,0,0,0,0,0), c(3, 2, N))
Y_M = Y

for (i in 1:N){
  Y[,1,i] = t(X1_2[i,])
  Y[,2,i] = t(X2_2[i,])
  Y_M[,1,i] = t(X1_2_M[i,])
  Y_M[,2,i] = t(X2_2_M[i,])
}

data = Y
data_M = Y_M

output_file = "./vcg_data/vcg_output_group2.RData"
data_with_init_with_MCMC_samples = finiteMixtureML(data,nc=1,true_nc=1,max_iter,run_id=1,output_file,vague_prior=3)

##### group 3 data
X1_3 = vcg_data[grp_3_id,3:5]
X2_3 = vcg_data[grp_3_id,6:8]

X1_3_M = vcg_data[grp_3_id,12:14]
X2_3_M = vcg_data[grp_3_id,15:17]

N = length(grp_3_id)
Y = array(c(0,0,0,0,0,0), c(3, 2, N))
Y_M = Y
  
  
for (i in 1:N){
  Y[,1,i] = t(X1_3[i,])
  Y[,2,i] = t(X2_3[i,])
  Y_M[,1,i] = t(X1_3_M[i,])
  Y_M[,2,i] = t(X2_3_M[i,])
}

data = Y
data_M = Y_M
output_file = "./vcg_data/vcg_output_group3.RData"
data_with_init_with_MCMC_samples = finiteMixtureML(data,nc=1,true_nc=1,max_iter,run_id=1,output_file,vague_prior=3)


##### group 4 data
X1_4 = vcg_data[grp_4_id,3:5]
X2_4 = vcg_data[grp_4_id,6:8]

X1_4_M = vcg_data[grp_4_id,12:14]
X2_4_M = vcg_data[grp_4_id,15:17]

N = length(grp_4_id)
Y = array(c(0,0,0,0,0,0), c(3, 2, N))
Y_M = Y

for (i in 1:N){
  Y[,1,i] = t(X1_4[i,])
  Y[,2,i] = t(X2_4[i,])
  Y_M[,1,i] = t(X1_4_M[i,])
  Y_M[,2,i] = t(X2_4_M[i,])
}

data = Y
data_M = Y_M
output_file = "./vcg_data/vcg_output_group4.RData"
data_with_init_with_MCMC_samples = finiteMixtureML(data,nc=1,true_nc=1,max_iter,run_id=1,output_file,vague_prior=3)




#########

### code below should give same output; use that in order to debug

#########
# set.seed(43185)
# print(date())
# 
# 
# N = dim(Y)[[3]]
# max_iter = 10
# MCMC_sample = vector("list", max_iter)
# data_with_init_with_MCMC_samples = NULL	
# data_with_init_with_MCMC_samples$data = data 
# nc = 1
# #### initialization ####
# #Rprof(interval = 0.02)
# ###curr_param = init_param_true(nc)
# init_param = init_param_from_MLE(data,nc)
# data_with_init_with_MCMC_samples$init_param = init_param
# #save(init_param_MLE=init_param,file=init_param_output_file)
# ########################
# print(date())
# ########################
# vague_prior = 1
# hyper = hyper_selection(nc,vague_prior,init_param)
# ########################
# 
# curr_param = init_param
# Z = curr_param$id_arr
# 
# #   prob_vec = rep(1,nc)
# #   prob_vec = prob_vec/nc
# #   for (i in 1:N){
# #     r = rmultinom(1,1,prob_vec)
# #     Z[i] = which(r==1)
# #   }
# 
# ### prior for mixture weights
# pi_vec = rep(1,nc)/nc
# #pi_vec = rdirichlet(1,hyper$dir_alpha) 
# cluster_id_cnt = sapply(1:nc, function(x) sum(curr_param$id_arr == x))
# pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)
# 
# d_tmp_val = matrix(rep(0,N*nc),ncol=nc)
# 
# for(iter in 1:max_iter){
#   ######################
#   if(iter%%10 == 0){
#     print(paste0("MCMC main iter = ",iter))
#   }
#   ##### gibbs step #####
#   #     t1=Sys.time()
#   #     clust_assign_prob_mat_0 = matrix(rep(0,N*nc),nrow=N)
#   #     for (i in 1:N){
#   #       clust_assign_prob_i = rep(0,nc)
#   #       for (cluster_id in 1:nc){
#   #         d_tmp_val[i,cluster_id] = dMatrixLangevin(data[,,i],curr_param$M[,,cluster_id],curr_param$D[,,cluster_id],curr_param$V[,,cluster_id])
#   #         
#   #         clust_assign_prob_i[cluster_id] = pi_vec[cluster_id]*d_tmp_val[i,cluster_id]
#   #       }
#   #       clust_assign_prob_mat_0[i,] = clust_assign_prob_i
#   #       #clust_assign_prob_i = clust_assign_prob_i/sum(clust_assign_prob_i)
#   #       Z[i] = sample(1:nc,1,FALSE,clust_assign_prob_i)
#   #       #print(clust_assign_prob_i)
#   #     }
#   #     Sys.time()-t1
#   #     
#   ### compact way to run the above commented code
#   t2=Sys.time()
#   clust_assign_prob_mat = matrix(rep(0,N*nc),nrow=N)
#   for (cluster_id in 1:nc){
#     d_tmp_val[,cluster_id] = apply(data,3,function(x) dMatrixLangevin(x,curr_param$M[,,cluster_id],curr_param$D[,,cluster_id],curr_param$V[,,cluster_id]))
#   }    
#   clust_assign_prob_mat = d_tmp_val*t(matrix(rep(pi_vec,N),ncol=N))
#   Z = apply(clust_assign_prob_mat,1,function(x) sample(1:nc,1,FALSE,x))
#   #Sys.time()-t2
#   
#   
#   cluster_assign_vec = Z
#   #write.table(d_tmp_val,file="d_tmp_val.txt")
#   #tmp_table = as.matrix(table(Z))
#   #cluster_id_cnt = cbind(as.integer(row.names(tmp_table)),tmp_table)
#   #cluster_id_cnt = rep(0,length(hyper$dir_alpha))
#   
#   #for (cluster_id in 1:nc){
#   #  cluster_id_cnt[cluster_id] = sum(cluster_assign_vec == cluster_id)
#   #}
#   ## as an alternative to previous 3 statements
#   cluster_id_cnt = sapply(1:nc, function(x) sum(cluster_assign_vec == x))
#   
#   #print(hyper$dir_alpha)
#   print(cluster_id_cnt)
#   pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)
#   
#   
#   ##############
# 
#   curr_param = cluster_param_update(Y, curr_param, cluster_assign_vec, cluster_id_vec=c(1), hyper)
# 
#   print(paste0("VCG one cluster modeling: MCMC iteration update done ",iter))
#   print(Sys.time()-t2)
#   MCMC_sample[[iter]]$curr_param = curr_param
#   MCMC_sample[[iter]]$pi_vec = pi_vec
#   
#   if(iter%%10==0){
#     data_with_init_with_MCMC_samples$MCMC_sample = MCMC_sample 
#     save(data_with_init_with_MCMC_samples, file = output_file)
#   }
# }
#   
# 
# 


##### Parameter Estimation for Part 1 Paper
analyse_MCMC_sample_for_vcg <-function(grp_id,burn_in){
  
  src_name = sprintf("./vcg_data/vcg_output_group%d.RData",grp_id)
  load(src_name)
  
  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample
  
  n_sam = length(MCMC_sample) - burn_in
  
  F_avg = (MCMC_sample[[1]]$curr_param$F)*0

  
  for (i in (burn_in+1):n_sam){
    
    curr_param = MCMC_sample[[i]]$curr_param
    
    F_avg = F_avg + curr_param$F
  }

  F_avg = F_avg/n_sam
  

  Downs_M_hat = matrix(c(0.770,0.623,0.136,0.606,-0.782,0.149),nrow = 3) 
  Downs_K_hat = matrix(c(6.51,0.36,0.36,13.40),nrow = 2)
}



