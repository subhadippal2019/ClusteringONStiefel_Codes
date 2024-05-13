#### run vcg for a single group
source("utility.R")
load_src_libs()

args = commandArgs(TRUE)

### k = as.numeric(args[1])
k1 = as.numeric(args[1])
k2 = as.numeric(args[2])

### find q and p from McFee systems

vcg_data = read.csv("./vcg_data/vcg.csv")
### for each of the four age groups

### grp_id = which(vcg_data$AgeSex==k)

grp_id1 = which(vcg_data$AgeSex==k1)
grp_id2 = which(vcg_data$AgeSex==k2)

grp_id = c(grp_id1,grp_id2)

max_iter = 5000

X1_M = vcg_data[grp_id,12:14]
X2_M = vcg_data[grp_id,15:17]

N = length(grp_id)

cat("total number of data in grp ",k1," and grp ",k2," is equal to ",N,"\n")

Y_M = array(c(0,0,0,0,0,0), c(3, 2, N))


for (i in 1:N){
 
  Y_M[,1,i] = t(X1_M[i,])
  Y_M[,2,i] = t(X2_M[i,])
}

data = Y_M

### output_file = sprintf("./vcg_data/vcg_output_group%d_McFee.RData",k)
output_file = sprintf("./vcg_data/vcg_output_group_%d_and_%d_McFee.RData",k1,k2)

data_with_init_with_MCMC_samples = finiteMixtureML(data,nc=1,true_nc=1,max_iter,run_id=1,output_file,vague_prior=3)

#### ========================

###analyse_MCMC_sample_for_vcg_one_grp <-function(grp_id,burn_in){

analyse_MCMC_sample_for_vcg_one_grp <-function(grp_id1,grp_id2,burn_in){
    
  #src_name = sprintf("./vcg_data/vcg_output_group%d_McFee.RData",grp_id)
  
  src_name = sprintf("./vcg_data/vcg_output_group_%d_and_%d_McFee.RData",grp_id1,grp_id2)
  
  load(src_name)

  MCMC_sample = data_with_init_with_MCMC_samples$MCMC_sample

  n_sam = length(MCMC_sample) - burn_in

  F_avg = (MCMC_sample[[1]]$curr_param$F)*0

  for (i in (burn_in+1):n_sam){

    curr_param = MCMC_sample[[i]]$curr_param

    F_avg = F_avg + curr_param$F
  }

  F_avg = F_avg/(n_sam-burn_in)

  return(F_avg)
  
  #Downs_M_hat = matrix(c(0.770,0.623,0.136,0.606,-0.782,0.149),nrow = 3)
  #Downs_K_hat = matrix(c(6.51,0.36,0.36,13.40),nrow = 2)

  #Downs_F_hat = Downs_M_hat%*%Downs_K_hat
  #Downs_F_hat
}


#### plot of 6 elements of matrix from group 1
grp_id = 1
src_name = sprintf("./vcg_data/vcg_output_group%d_McFee.RData",grp_id)
load(src_name)
i=1
for(i in 1:1){
  T=data_with_init_with_MCMC_samples$MCMC_sample[[i]]$curr_param
  M=matrix(T$M,nrow=3)
  D=matrix(T$D,nrow=2)
  V=matrix(T$V,nrow=2)
  
  
}