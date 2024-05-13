

#R CMD BATCH --no-save --no-restore '--args ./ML_dataset_4_2.Rdata 2' ./runAllCode.R ./run2.log

rm(list=ls(all=TRUE))
args <- commandArgs(TRUE)
set.seed(72357)
source("utility.R")
load_src_libs()

#L = generate_simulated_data_ML(N=300)
load(args[1])
n_clust_to_run = as.numeric(args[2])

for(i in n_clust_to_run){
  
  nClust = i

  
  ulist_str1 = unlist(strsplit(args[1],"[/]"))
  L1=length(ulist_str1) 
  str1 = ulist_str1[L1]
  true_nClust = unlist(strsplit(str1,"[_.]"))[3]
  data_id = unlist(strsplit(str1,"[_.]"))[4]

  ##data,nc,max_iter=2,run_id=1 
  output_file = sprintf("new_output_nclust_%d_data_id_%s_true_nClust_%s.RData",nClust,data_id,true_nClust)
  data_with_MCMC_sample = finiteMixtureML(L$data,nClust,true_nClust,max_iter=20,run_id=as.numeric(data_id),output_file)
  
  
  #data_with_MCMC_sample = NULL
  #data_with_MCMC_sample$L = L
  #data_with_MCMC_sample$MCMC_sample = MCMC_sample
  #save(data_with_MCMC_sample,file=out_file_name)
  #load("MCMC_sample.Rdata")
  #analyse_MCMC_sample(MCMC_sample)
  #plot_mcmc_sample_D(MCMC_sample,L)
  
}  
  
