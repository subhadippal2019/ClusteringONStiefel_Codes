#' parametric mixture modeling of Matrix Langevin with fixed clusters
#' @export
#' @param data input-data coming from Stiefel array of 2D matrices data[,,i] is in V_{3,2}
#' @param nc numer of clusters to fit the model
#' @param max_iter total number of MCMC ierations
#' @useDynLib ClusterOnStiefel
finiteMixtureML <- function(data,nc,max_iter=2,run_id=1){

  #library(gtools)

  set.seed(43185)

  N = dim(data)[[3]]

  vague_prior = 1
  hyper = init_run(nc,vague_prior)
  ########################

  MCMC_sample = vector("list", max_iter)

  MCMC_output_file = sprintf("MCMC_sample_%d_%d.RData",nc,run_id)
  init_param_output_file = sprintf("init_param_MLE_%d_%d.RData",nc,run_id)
  ########### hyper parameters #######
  ### need to select empirical prior
#   hyper=NULL
#   hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
#   hyper$H = (matrix(c(0,0,0,0),ncol = 2))
#   hyper$alpha = 1.0
#   hyper$beta = 0.0
#   hyper$dir_alpha = rep(0.0,nc)
#   for (j in 1:nc)
#     hyper$dir_alpha[j] = 1.0
#
  ####################################

  #### initialization ####

  ###curr_param = init_param_true(nc)
  init_param = init_param_from_MLE(data,nc)
  save(init_param_MLE=init_param,file=init_param_output_file)
  ########################

  curr_param = init_param
  Z = curr_param$id_arr

#   prob_vec = rep(1,nc)
#   prob_vec = prob_vec/nc
#   for (i in 1:N){
#     r = rmultinom(1,1,prob_vec)
#     Z[i] = which(r==1)
#   }

  ### prior for mixture weights
  pi_vec = rep(1,nc)/nc
  #pi_vec = rdirichlet(1,hyper$dir_alpha)
  cluster_id_cnt = sapply(1:nc, function(x) sum(curr_param$id_arr == x))
  pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)

  d_tmp_val = matrix(rep(0,N*nc),ncol=nc)

  for(iter in 1:max_iter){
    ######################
    if(iter%%10 == 0){
      print(paste0("MCMC main iter = ",iter))
    }
    ##### gibbs step #####
    for (i in 1:N){
      clust_assign_prob_i = rep(0,nc)
      for (cluster_id in 1:nc){
        d_tmp_val[i,cluster_id] = dMatrixLangevin(data[,,i],curr_param$M[,,cluster_id],curr_param$D[,,cluster_id],curr_param$V[,,cluster_id])

        clust_assign_prob_i[cluster_id] = pi_vec[cluster_id]*d_tmp_val[i,cluster_id]
      }
      #clust_assign_prob_i = clust_assign_prob_i/sum(clust_assign_prob_i)
      Z[i] = sample(1:nc,1,FALSE,clust_assign_prob_i)
      #print(clust_assign_prob_i)
    }
    cluster_assign_vec = Z
    #write.table(d_tmp_val,file="d_tmp_val.txt")
    #tmp_table = as.matrix(table(Z))
    #cluster_id_cnt = cbind(as.integer(row.names(tmp_table)),tmp_table)
    #cluster_id_cnt = rep(0,length(hyper$dir_alpha))

    #for (cluster_id in 1:nc){
    #  cluster_id_cnt[cluster_id] = sum(cluster_assign_vec == cluster_id)
    #}
    ## as an alternative to previous 3 statements
    cluster_id_cnt = sapply(1:nc, function(x) sum(cluster_assign_vec == x))

    #print(hyper$dir_alpha)
    print(cluster_id_cnt)
    pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)

    #### for all the clusters ####
    #if(iter == 1){
    #  load("ML_dataset.Rdata")
    #  curr_param = L$curr_param
    #  cluster_assign_vec = L$clust
    #}
    curr_param = cluster_param_update(data, curr_param, cluster_assign_vec, cluster_id_vec=1:nc,hyper)

    ### need to update with updated cluster_assign_vec
    ##################################

    print(paste0("Mixture modeling: MCMC iteration update ",iter))

    MCMC_sample[[iter]]$curr_param = curr_param
    MCMC_sample[[iter]]$pi_vec = pi_vec

    if(iter%%100==0){
      save("MCMC_sample", file = MCMC_output_file)
    }
  }

  #G = NULL
  #G$curr_param = curr_param
  #G$clust = cluster_assign_vec
  #print(G$curr_param$D[,,1])

  #return(G)
  #save(D1_1,D2_1,D1_2,D2_2,D1_3,D2_3,file="D1_D2.Rdata")
  save("MCMC_sample", file = MCMC_output_file)

  return(MCMC_sample)
}


dMatrixLangevin <- function(X,M,D,V){

  #dyn.load("./C_code/functions_hyper_2by2_R.so")

  F = M %*% D %*% V
  D = diag(D)
  eigenValues=D^2/4;
  dRet = 0.0
  ### normalizing constant
  hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,dRet)[[3]]

  ### compute density value
  ML_density =exp(sum(diag(t(F)%*%X))) / hyper0F1_val

  return(ML_density)

}


init_run <-function(nc,vague_prior){
  source("utility.R")
  load_src_libs()

  hyper=NULL

  if(vague_prior == 1){
    hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
    hyper$H = (matrix(c(0,0,0,0),ncol = 2))
    hyper$alpha = 1.0
    hyper$beta = 0.15 ###### note that, min 5 data made a cluster 5*0.01 = 0.05,
                    ###### for N==1 if S>0.99 then distribution of d1 hv a heavy tail

    hyper$dir_alpha = rep(0.0,nc)
    for (j in 1:nc)
      hyper$dir_alpha[j] = 1.0
  } ## emperical


  hyper$debug = 0

  return(hyper)

}
