source("utility.R")
load_src_libs()
source("analyse_MCMC_sample.R")
result_NEC = NULL

for (case in c(1,2,3)){
for(k in c(500,1000,1500)){
for (i in 3:6){
  fileName = sprintf("../RealData/NEC/output/output_nclust_%d_realNEC.RData",i)
  tmp1 = strsplit(fileName,"/")[[1]][5]
  str1 = sprintf("%s_%s_%s",strsplit(tmp1,"_")[[1]][2],strsplit(tmp1,"_")[[1]][3],strsplit(tmp1,"_")[[1]][4])
  result_NEC$sample = str1
  load(fileName)
  burn_in = k
  
  result_NEC$nclust = i
  if((case == 1)){
    res = calculate_DIC_5(data_with_MCMC_sample,burn_in)
    result_NEC$val = res$val
    result_NEC$pd = res$pd
    file_name = sprintf("DIC_5_NEC_%d_2000_mcmc_samples.txt",burn_in)
  }
  
  if((case == 2)){
    res = calculate_DIC_5_with_penalty(data_with_MCMC_sample,burn_in)
    result_NEC$val = res$val
    result_NEC$pd = res$pd
    file_name = sprintf("DIC_NEC_%d_with_penalty_2000_mcmc_samples.txt",burn_in)
  }
  
  if((case == 3)){
    res = calculate_DIC(data_with_MCMC_sample$MCMC_sample,data_with_MCMC_sample$data,burn_in)
    result_NEC$val = res
    #result_NEC$pd = res$pd
    file_name = sprintf("DIC_NEC_%d_original_2000_mcmc_samples.txt",burn_in)
  }
  
  write.table(result_NEC,sep="\t",quote=F,row.name=F,col.names=F,append=T,file=file_name)
}
}
}

