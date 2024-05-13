
library(gsl)
library(Rcpp)
sourceCpp('C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/package_code/src/functions_hyper_2by2_CPP.cpp')

#hyper_2by2_R(1.5,c(3,1))

P_Data_given_param_ML<-function(Data, M,D,V, log=T){
 n_data=dim(Data)[3]
 p=dim(Data)[2]
 D= matrix(D,ncol=p);V= matrix( V, ncol=p);M= matrix( M, ncol=p);
F=M%*%D%*%V
 X_bar=apply(Data,c(1,2), mean)
  log_numerator=sum( diag(    (t(F)%*%X_bar )    ))
  denominator= hyper_2by2_R(1.5,diag(D))
  log_P= n_data*( log_numerator-  log(denominator) )
  if(log){
    val=log_P
  }
  if(!log){
    val=exp(log_P)
  }
  return (val)
  
}


# 
# P_Data_given_param_ML<-function(Data, M,D,V, log=T){
#   n_data=dim(Data)[3]
#   p=dim(Data)[2]
#   D= matrix(D,ncol=p);V= matrix( V, ncol=p);M= matrix( M, ncol=p);
#   
#   X_bar=apply(Data,c(1,2), mean)
#   log_numerator=sum( diag(    (V%*%D%*%t(M)%*%X_bar )    ))
#   denominator= hyper_2by2_R(1.5,diag(D))
#   log_P= n_data*( log_numerator-  log(denominator) )
#   if(log){
#     val=log_P
#   }
#   if(!log){
#     val=exp(log_P)
#   }
#   return (val)
#   
# }

P_Data_given_param_ML_from_MCMC<-function(Data,MCMC,burnIN=0, log=T ){
  
  MC_length=length(MCMC)
  start=burnIN+1;
  val=P_Data_given_param_ML(Data,  MCMC[[start]]$curr_param$M,MCMC[[start]]$curr_param$D,MCMC[[start]]$curr_param$V,log)
  for (i  in (start+1):MC_length){
  temp=P_Data_given_param_ML(Data,  MCMC[[i]]$curr_param$M,MCMC[[i]]$curr_param$D,MCMC[[i]]$curr_param$V)
  val=c(val,temp)
  }
  
  return(val)
}



Reciprocal_BayesFactorConstant<-function(val){
  
  #The intendent normalizing constant is 1/(exp(Log_multi)*MeanNormalizedVal)
  
  multiplier=max(-val)
  val_normalized=(-val)-multiplier
  MeanNormalizedVal=mean(exp( val_normalized))
  lst=list(reciproCalMean=MeanNormalizedVal, Log_multi=multiplier)
  return(lst)
  
}



log_Reciprocal_BayesFactorConstant<-function(val){
  
  
  temp_val=-val
  #The intendent normalizing constant is 1/(exp(Log_multi)*MeanNormalizedVal)
  Max=max(temp_val)
  log_Z=log(mean(exp(temp_val-Max)))+ Max
  return( log_Z)
  
}

  
Calculate_ByesFactor<-function(Gr1,Gr2,burnIN=3000){
#  browser()
  FileLocation="C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/"
  FileName1= paste0("vcg_output_group",Gr1,"_McFee.RData")
  FileName2= paste0("vcg_output_group",Gr2,"_McFee.RData")
  FileName3= paste0("vcg_output_group_",Gr1,"_and_",Gr2,"_McFee.RData")
  load(paste0(FileLocation,FileName1))
  val_Gr1=P_Data_given_param_ML_from_MCMC(data_with_init_with_MCMC_samples$data,data_with_init_with_MCMC_samples$MCMC_sample,burnIN)
  
  load(paste0(FileLocation,FileName2))
  val_Gr2=P_Data_given_param_ML_from_MCMC(data_with_init_with_MCMC_samples$data,data_with_init_with_MCMC_samples$MCMC_sample,burnIN)
  
  
  load(paste0(FileLocation,FileName3))
  val_Gr1_Gr2=P_Data_given_param_ML_from_MCMC(data_with_init_with_MCMC_samples$data,data_with_init_with_MCMC_samples$MCMC_sample,burnIN)
  
  
  #Ratio_1=Reciprocal_BayesFactorConstant(val_Gr1)
  
#  Ratio_2=Reciprocal_BayesFactorConstant(val_Gr2)
  
 # Ratio_1_2=Reciprocal_BayesFactorConstant(val_Gr1_Gr2)
  
  #fractionRatio=(Ratio_1$reciproCalMean*Ratio_2$reciproCalMean/Ratio_1_2$reciproCalMean)
  #fractionLogMult=Ratio_1$Log_multi+Ratio_2$Log_multi-Ratio_1_2$Log_multi
  #log_BayesFactor= log(fractionRatio)+fractionLogMult
  
  Log_Ratio_1=-log_Reciprocal_BayesFactorConstant(val_Gr1)
  Log_Ratio_2=-log_Reciprocal_BayesFactorConstant(val_Gr2)
  Log_Ratio_3=-log_Reciprocal_BayesFactorConstant(val_Gr1_Gr2)
  log_BayesFactor=Log_Ratio_3-(Log_Ratio_1+Log_Ratio_2)
  
  return(log_BayesFactor)
  
}


bfl1_3=Calculate_ByesFactor(1,3,4000)
bfl2_4=Calculate_ByesFactor(2,4,4000)
# 
pdf(file='Plot_Logliklihood_Group1_BAyesFactor.pdf', height=5 , width=6)
plot(val_1,main='Loglikelihood evaluated at different MCMC sample for Group 1',ylab='Likelihood', xlab='Index for Markov chain sample')
dev.off()

  
  pdf(file='Plot_Logliklihood_Group3_BAyesFactor.pdf', height=5 , width=6)
  plot(val_3,main='Loglikelihood evaluated at different MCMC sample for Group 3',ylab='Likelihood', xlab='Index for Markov chain sample')
  dev.off()
  
  
  pdf(file='Plot_Logliklihood_Group1_3_BAyesFactor.pdf', height=5 , width=6)
  plot(val_1_3,main='Loglikelihood evaluated at different MCMC sample for Group 1 and Group 3 combinde',ylab='Likelihood', xlab='Index for Markov chain sample')
  dev.off()
# start=10
# MM= MCMC[[start]]$curr_param$M
# VV= MCMC[[start]]$curr_param$V
# DD= MCMC[[start]]$curr_param$D
# FF= MCMC[[start]]$curr_param$F
# matrix(FF,ncol=2)-matrix(MM,ncol=2)%*%matrix(DD,ncol=2)%*%(matrix(V,ncol=2))
# 
# 
# load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group1_McFee.RData")
# val_1=P_Data_given_param_ML_from_MCMC(data_with_init_with_MCMC_samples$data,data_with_init_with_MCMC_samples$MCMC_sample)
# 
# load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group_1_and_3_McFee.RData")
# val_1_3=P_Data_given_param_ML_from_MCMC(data_with_init_with_MCMC_samples$data,data_with_init_with_MCMC_samples$MCMC_sample)
# 
# 
# load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group1_McFee.RData")
# val_1=P_Data_given_param_ML_from_MCMC(data_with_init_with_MCMC_samples$data,data_with_init_with_MCMC_samples$MCMC_sample)
# 
# load("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/vcg_data/vcg_output_group3_McFee.RData")
# val_3=P_Data_given_param_ML_from_MCMC(data_with_init_with_MCMC_samples$data,data_with_init_with_MCMC_samples$MCMC_sample)
