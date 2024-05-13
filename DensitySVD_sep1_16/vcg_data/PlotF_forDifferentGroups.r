


extract_F_frm_MCMC<-function(data_with_init_with_MCMC_samples, burnIN=0,grp=1){
  McLen=length(data_with_init_with_MCMC_samples$MCMC_sample)
  F=array(0,dim=c(3,2,McLen-burnIN))
  for(mcINdex in 1:(McLen-burnIN)){
    F[,,mcINdex]=data_with_init_with_MCMC_samples$MCMC_sample[[mcINdex+burnIN]]$curr_param$F[,,1]
    
    
  }
File=  'C:\\Users\\subha\\Dropbox\\projects\\ClusteringDTIonStiefel\\DensitySVD_sep1_16\\vcg_data\\F_summary\\EstimatedDensityForF_Group_'
write.csv(cbind(MeanF=lst$MeanF,SdF=lst$SdF),file=paste(File,grp,'_burnIN2000','.csv'))
  pdf(paste(File,grp,'_burnIN2000','.pdf'),width = 9, height=6)
  par(mfrow=c(2,3))
  for(i in 1:3){
    for(j in 1:2){
      plot(density(F[i,j,]),col='blue',main=paste0('Density Plot for (',i,',',j,') th component of F' ), ylab='Estimated Density', xlab='Magnitude')
    }
  }
  dev.off()
  
 listF= list(F=F,MeanF=apply(F,c(1,2),mean),SdF=apply(F, c(1,2),sd))
  
  return(listF)
}

load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group1_McFee.RData")
lst=extract_F_frm_MCMC(data_with_init_with_MCMC_samples, 2000,1)

load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group2_McFee.RData")
lst=extract_F_frm_MCMC(data_with_init_with_MCMC_samples, 2000,2)

load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group3_McFee.RData")
lst=extract_F_frm_MCMC(data_with_init_with_MCMC_samples, 2000,3)

load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group4_McFee.RData")
lst=extract_F_frm_MCMC(data_with_init_with_MCMC_samples, 2000,4)
