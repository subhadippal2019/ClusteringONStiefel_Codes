
A=read.csv("DIC_table_for_nclust_3.csv")
dic_id_3 = apply(A,1,which.min)
dic_cnt_3 = sum(dic_id_3==3)

B=read.csv("DIC_table_for_nclust_4.csv")
dic_id_4 = apply(B,1,which.min)
dic_cnt_4 = sum(dic_id_4==4)

DIC_VERSION = 5 ## 4,5,6,7,8

COL_ID = (DIC_VERSION-3)*2 


sim3 = read.csv("sim_true_3_dic_4_5_6_7_8.txt",sep='\t',header=F)
mat3 = matrix(rep(0,200),nrow=50)
for (dId in 1:50){
  k1 = dId
  k2 = dId + 50
  k3 = dId + 100 
  k4 = dId + 150
  mat3[dId,] = c(sim3[k1,COL_ID],sim3[k2,COL_ID],sim3[k3,COL_ID],sim3[k4,COL_ID])
  
}
id3 = apply(mat3,1,which.min)


##########################

sim4 = read.csv("sim_true_4_dic_4_5_6_7_8.txt",sep='\t',header=F)
mat4 = matrix(rep(0,250),nrow=50)
for (dId in 1:50){
  k1 = dId
  k2 = dId + 50
  k3 = dId + 100 
  k4 = dId + 150
  k5 = dId + 200
  mat4[dId,] = c(sim4[k1,COL_ID],sim4[k2,COL_ID],sim4[k3,COL_ID],sim4[k4,COL_ID],sim4[k5,COL_ID])
  
}

id4 = apply(mat4,1,which.min)

cat('with DIC for 3 : cases = ',dic_cnt_3,'\n')
cat('with DIC for 4 : cases = ',dic_cnt_4,'\n')


cat('with DIC_5 for 3 : cases = ',sum(id3==2),'\n')
cat('with DIC_5 for 4 : cases = ',sum(id4==3),'\n')


##########################
#### summary for F-F_hat
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


calc_MSE_F__F_hat <- function(){
  library(combinat)
  data_set_len = 50
  MSE_perm_matrix = NULL
  MSE_matrix = matrix(rep(0,2*data_set_len),ncol=data_set_len)
  Perm_matrix = matrix(rep(0,2*data_set_len),ncol=data_set_len)
  Rel_change_matrix = matrix(rep(0,2*data_set_len),ncol=data_set_len)
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
      Rel_change_matrix[t_clust-2,i] = ret_MSE$val/total_cluster_norm(F_t^2)
    }
  }
  MSE_perm_matrix$MSE_matrix = MSE_matrix
  MSE_perm_matrix$Perm_matrix = Perm_matrix
  MSE_perm_matrix$MSE_relative_change_matrix = Rel_change_matrix
  return(MSE_perm_matrix)
}

total_cluster_norm <-function(F_t){
  nClust = dim(F_t)[3]
  S = 0.0
  for(i in 1:nClust){
    S = S+norm(F_t[,,i],type="F")
  }
  S = S/(nClust)
  
  return(S)
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
      #S = S+norm(F_t[,,i] - F_h_perm[,,i])
      
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

MSE_perm_matrix = calc_MSE_F__F_hat()

pdf("F_minus_F_hat_rel_change.pdf",height=5,width=7)
par(mfrow=c(1,2))
plot(c(MSE_perm_matrix$MSE_relative_change_matrix[1,])*100,
      xlab="dataset index", ylab=expression('d'[MSE]*' as a percentage of absolute sum of elements in F'),pch=19,ylim=c(0,5),xaxt='n',yaxt='n')
abline(h=4,col='blue',lty=2)
axis(side = 1, at = c(1,10,20,30,40,50),labels=c(1,10,20,30,40,50))
axis(side = 2, at = c(1,2,3,4,5),labels=c('1%','2%','3%','4%','5%'))

plot(c(MSE_perm_matrix$MSE_relative_change_matrix[2,])*100,
     xlab="dataset index", ylab = "",pch=19,ylim = c(0,5),xaxt='n',yaxt='n')
abline(h=4,col='blue',lty=2)
axis(side = 2, at = c(1,2,3,4,5),labels=c('1%','2%','3%','4%','5%'))

#plot(c(MSE_perm_matrix$MSE_relative_change_matrix[2,])*100,
#     xlab="dataset index", ylab=expression("Relative change between true F and estimated F (in percentage)"),pch=19,cex=0.5,xaxt='n')

axis(side = 1, at = c(1,10,20,30,40,50),labels=c(1,10,20,30,40,50))
axis(side = 2, at = c(1,2,3,4,5),labels=c('1%','2%','3%','4%','5%'))
#abline(v=50.5,col = 'blue')
dev.off()

