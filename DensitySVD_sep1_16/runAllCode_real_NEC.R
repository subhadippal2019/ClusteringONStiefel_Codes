

#R CMD BATCH --no-save --no-restore '--args ./real_DTI_data.Rdata 6' ./runAllCode_real.R ./run_real_14206_100.log
#R CMD BATCH --no-save --no-restore '--args ./nec_175_data.Rdata 6' ./runAllCode_real_NEC.R ./run_real_nec_175.log

#R CMD BATCH --no-save --no-restore '--args ./nec_175_data_july27.RData 6' ./runAllCode_real_NEC.R ./run_real_nec_175_july28.log

rm(list=ls(all=TRUE))

args <- commandArgs(TRUE)
source("finiteMixtureML.R")
set.seed(77357)

#L = generate_simulated_data_ML(N=300)
load(args[1])
n_clust_to_run = as.numeric(args[2])

for(i in 3:n_clust_to_run){
  
  nClust = i

  #true_nClust = strsplit(args[[1]],"[_.]")[[1]][4]
  #data_id = 12345
 
  ##data,nc,max_iter=2,run_id=1 
  #finiteMixtureML <- function(data,nc,true_nc,max_iter=2,run_id=1,output_file){
  
  #out_file_name = sprintf("output_nclust_%d_realDTI_%d.RData",nClust,data_id)
  #out_file_name = sprintf("../RealData/NEO/output/output_nclust_%d_realNEO_%d.RData",nClust,data_id)
  out_file_name = sprintf("../RealData/NEC/output/july27_output_nclust_%d_realNEC.RData",nClust)
  
  #data_with_MCMC_sample = finiteMixtureML(data=neo_192,nClust,nClust,max_iter=1000,run_id=as.numeric(data_id),output_file = out_file_name)
  data_with_MCMC_sample = finiteMixtureML(data=nec_175,nClust,nClust,max_iter=2000,run_id=as.numeric(data_id),output_file = out_file_name)
  
  
  #data_with_MCMC_sample = NULL
  #L = neo_192
  #data_with_MCMC_sample$L = L
  #data_with_MCMC_sample$MCMC_sample = MCMC_sample
  save(data_with_MCMC_sample,file=out_file_name)
  #load("MCMC_sample.Rdata")
  #analyse_MCMC_sample(MCMC_sample)
  #plot_mcmc_sample_D(MCMC_sample,L)
  
}  
  
