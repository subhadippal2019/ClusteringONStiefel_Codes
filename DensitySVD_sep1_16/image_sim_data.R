#### simulation data generation with image information
library(spatstat)

nRow = 40; nCol= 40
mat_indx = matrix(rep(0,nRow*nCol),ncol=nCol)

load("~/Desktop/DTI_project/simulated_data/ML_dataset_3_1.RData")
T=table(L$clust)
img_r = 40
img_c = 40

assign_mat = matrix(rep(0,nRow*nCol),ncol=nCol)

par(mfrow=c(2,1))

init_val = c(3,18,26,35)

for(i in 1:length(T)){
  mat_indx[(init_val[i]*img_r+1):(init_val[i]*img_r+T[i])] = i
  assign_mat[(init_val[i]*img_r+1):(init_val[i]*img_r+T[i])] = which(L$clust==i)  
}
plot(im(mat_indx))


load("~/Desktop/DTI_project/simulated_data/results_nclust_3/outputs/output_nclust_3_data_id_1_true_nClust_3.RData")
est_id_arr = data_with_init_with_MCMC_samples$MCMC_sample[[1000]]$curr_param$id_arr
mat_indx_prime = mat_indx

mat_indx_prime[which(assign_mat>0)] = est_id_arr[assign_mat[which(assign_mat>0)]]
plot(im(mat_indx_prime))
