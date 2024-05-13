
norm_vec <- function(x){
  return(sqrt(sum(x*x)))
}

computeFA <-function(x){

  if(length(x) != 3){
    stop("Error !!")
  }
  
  LocateInvalidity=which(x<=0)
  
  if(length(LocateInvalidity) > 0){
    if(length(LocateInvalidity)==length(x)){
      FA=0
    }
  
    if(length(LocateInvalidity)<length(x)){
      x[LocateInvalidity]=min(x[which(x>0)])*0.0000001
      #print(x)
      FA=sqrt(length(x))*sd(x)/norm_vec(x)
    }
  }
  if(length(LocateInvalidity) == 0){
    FA = sqrt(length(x))*sd(x)/norm_vec(x)
  }
  
  return (FA)
}

load_data_from_nii_file = function(){

  #imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_L1.nii")
  #imgL2=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_L2.nii")
  #imgL3=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_L3.nii")
  
  
  imgL1=readNIfTI("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/RealData/PNC/600085654611/data_L1.nii.gz")
  imgL2=readNIfTI("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/RealData/PNC/600085654611/data_L2.nii.gz")
  imgL3=readNIfTI("C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/RealData/PNC/600085654611/data_L3.nii.gz")
  
  
  L1 = imgL1@.Data
  L2 = imgL2@.Data
  L3 = imgL3@.Data
  
  Lambda=array(rep(0.0,prod(dim(L1))*3),c(dim(L1)[1],dim(L1)[2],dim(L1)[3],3))
  Lambda[,,,1]=L1
  Lambda[,,,2]=L2
  Lambda[,,,3]=L3
  
  Lambda1=array(Lambda, c(prod(dim(L3)),3))
  check=rep(0.0,prod(dim(L3)))
  #for (i in 1:prod(dim(L3))){
  #  check[i]=1-(Lambda1[i,1]>= Lambda1[i,2])*(Lambda1[i,2]>=Lambda1[i,3])
  #}

  #check = apply(Lambda1,1,function(x){1-((x[1]>=x[2])*(x[2]>=x[3]))})  

  print(paste0("before FA compute!!"))
  
  FA=rep(0.0,prod(dim(L3)))
  #for (i in 1:prod(dim(L3))){
  #  FA[i] = computeFA(Lambda1[i,])
  #}
  FA = apply(Lambda1,1,computeFA)
  
  save(FA,file='C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/RealData/PNC/600085654611/FA_new_PNC.RData')

#   Lambda_corrected = array(rep(0.0,3*prod(dim(L1))),c(3,prod(dim(L3))))
#   
#   for (i in 1:prod(dim(L3))){
#     Lambda_corrected[,i]=Lambda_correction(Lambda1[i,]);
#   }
#   save(Lambda_corrected,file='./Corrected_lambda_new.RData')
#   
  
}

find_validVoxel <- function(cutoff){
  load('./FA_new_PNC.RData')
  validVoxel=which(FA>cutoff)
  return(validVoxel)
}

#### this actually generates the real data once FA cut-off value is given
get_valid_vectors <-function(cutoff){
  library(oro.nifti)
  #imgV1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V1.nii")
  #imgV2=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V2.nii")
  #imgV3=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V3.nii")
  
  imgV1=readNIfTI("C:\\Users\\subha\\Dropbox\\projects\\ClusteringDTIonStiefel\\RealData\\PNC\\600085654611\\data_V1.nii.gz")
  imgV2=readNIfTI("C:\\Users\\subha\\Dropbox\\projects\\ClusteringDTIonStiefel\\RealData\\PNC\\600085654611\\data_V2.nii.gz")
  
  
  V1 = imgV1@.Data
  V2 = imgV2@.Data
  
  V1_1=array(V1, c(prod(dim(V1)[1],dim(V1)[2],dim(V1)[3]),3))
  V2_1=array(V2, c(prod(dim(V2)[1],dim(V2)[2],dim(V2)[3]),3))
  
  selected_indx = find_validVoxel(cutoff)
  L1 = length(selected_indx)
  
  print(paste0("length of selected data = ",L1))
  
  real_DTI_data = array(rep(0.0,3*2*L1),c(3,2,L1))
  
  for(i in 1:L1){
    k = selected_indx[i]
    real_DTI_data[,,i] = cbind(t(V1_1[k,]),t(V2_1[k,]))
  }
  
  fileName=sprintf('C:\\Users\\subha\\Dropbox\\projects\\ClusteringDTIonStiefel\\RealData\\PNC\\600085654611\\real_DTI_data_%0.2f_PNC_large.RData',cutoff)
  save(real_DTI_data,file=fileName)

}


##### CREATING MASK

map_back_to_brain <-function(cutoff=0.31,data_from_result,extn){
  #### data_from_result : RData
  
  library(oro.nifti)
  #imgV1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V1.nii")
  #imgV2=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V2.nii")
  #imgV3=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V3.nii")
  
  #V1 = imgV1@.Data
  #V2 = imgV2@.Data
  
  #V1_1=array(V1, c(prod(dim(V1)[1],dim(V1)[2],dim(V1)[3]),3))
  #V2_1=array(V2, c(prod(dim(V2)[1],dim(V2)[2],dim(V2)[3]),3))
  
  #imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_L1.nii")
  imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/PNC/600085654611/data_L1.nii.gz")
  
  #imgL1=readNIfTI("~/Desktop/data_L1.nii")
  L1 = imgL1@.Data
  L1_1=array(L1, c(prod(dim(L1)[1],dim(L1)[2],dim(L1)[3]),1))
  L1_1 = L1_1*0
  
  #load('~/Desktop/results_real_6_7_8_1000_beta/outputs/output_nclust_6_real_data_id_12345.RData')
  load(data_from_result)
  selected_indx = find_validVoxel(cutoff)
  mcmc_len = length(data_with_init_with_MCMC_samples$MCMC_sample)
  #mcmc_len=500
  print(mcmc_len)

  L1_1[selected_indx] = data_with_init_with_MCMC_samples$MCMC_sample[[mcmc_len]]$curr_param$id_arr
  #L1_1[selected_indx] = data_with_init_with_MCMC_samples$init_param$id_arr
  print(paste0("length of selected data = ",length(L1_1)))
  
  L1_2=array(L1_1, c(dim(L1)[1],dim(L1)[2],dim(L1)[3]))
  
  imgL1@.Data = L1_2
  
  fileName=sprintf('~/Desktop/gen_L1_%s_%s',cutoff,extn)
  writeNIfTI(imgL1,fileName)
  
}

map_back_to_brain_binary <-function(cutoff=0.31,cluster_id,data_from_result,extn){
  
  library(oro.nifti)
  #imgV1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V1.nii")
  #imgV2=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V2.nii")
  #imgV3=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V3.nii")
  
  #V1 = imgV1@.Data
  #V2 = imgV2@.Data
  
  #V1_1=array(V1, c(prod(dim(V1)[1],dim(V1)[2],dim(V1)[3]),3))
  #V2_1=array(V2, c(prod(dim(V2)[1],dim(V2)[2],dim(V2)[3]),3))
  
  #imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_L1.nii")
  imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/PNC/600085654611/data_L1.nii.gz")
  #imgL1=readNIfTI("~/Desktop/data_L1.nii")
  L1 = imgL1@.Data
  L1_1=array(L1, c(prod(dim(L1)[1],dim(L1)[2],dim(L1)[3]),1))
  L1_1 = L1_1*0
  
  #load('~/Desktop/results_real_6_7_8_1000_beta/outputs/output_nclust_6_real_data_id_12345.RData')
  load(data_from_result)
  cat("length(selected_indx) = ",length(selected_indx),"/n")
  mcmc_len = length(data_with_init_with_MCMC_samples$MCMC_sample)
  print(mcmc_len)
  mcmc_len=210
  #L1_1[selected_indx] = data_with_init_with_MCMC_samples$MCMC_sample[[mcmc_len]]$curr_param$id_arr
  L1_1[selected_indx] = cluster_id*(data_with_init_with_MCMC_samples$init_param$id_arr == cluster_id)
  print(paste0("length of selected data = ",length(L1_1)))
  
  L1_2=array(L1_1, c(dim(L1)[1],dim(L1)[2],dim(L1)[3]))
  
  imgL1@.Data = L1_2
  
  fileName=sprintf('~/Desktop/new_gen_L1_%s_%s_%s',cutoff,extn,cluster_id)
  writeNIfTI(imgL1,fileName)
  
}


map_back_to_brain_binary_all <-function(cutoff=0.31,soft_cluster_for_data){
  ### 
  
  library(oro.nifti)
  #imgV1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V1.nii")
  #imgV2=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V2.nii")
  #imgV3=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V3.nii")
  
  #V1 = imgV1@.Data
  #V2 = imgV2@.Data
  
  #V1_1=array(V1, c(prod(dim(V1)[1],dim(V1)[2],dim(V1)[3]),3))
  #V2_1=array(V2, c(prod(dim(V2)[1],dim(V2)[2],dim(V2)[3]),3))
  
  #imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_L1.nii")
  #imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/PNC/600085654611/data_L1.nii.gz")
  #imgL1=readNIfTI("~/Desktop/data_L1.nii")
  
  L1 = imgL1@.Data
  L1_1=array(L1, c(prod(dim(L1)[1],dim(L1)[2],dim(L1)[3]),1))
  L1_1 = L1_1*0
  
  #load('~/Desktop/results_real_6_7_8_1000_beta/outputs/output_nclust_6_real_data_id_12345.RData')
  #load(data_from_result)
  selected_indx = find_validVoxel(cutoff)
  #mcmc_len = length(data_with_init_with_MCMC_samples$MCMC_sample)
  #print(mcmc_len)
  #mcmc_len=210
  #L1_1[selected_indx] = data_with_init_with_MCMC_samples$MCMC_sample[[mcmc_len]]$curr_param$id_arr
  C = dim(soft_cluster_for_data)[2]
  for (cluster_id in 1:C){
    L1_1[selected_indx] = soft_cluster_for_data[,cluster_id]
    print(paste0("length of selected data = ",length(L1_1)))
    
    L1_2=array(L1_1, c(dim(L1)[1],dim(L1)[2],dim(L1)[3]))
    
    imgL1@.Data = L1_2
    
    fileName=sprintf('~/Desktop/new_gen_L1_%s_%d',cutoff,cluster_id)
    writeNIfTI(imgL1,fileName)
  }
  
}

map_back_to_brain_binary_all_1 <-function(cutoff=0.31,soft_cluster_for_data){
  
  library(oro.nifti)
  #imgV1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V1.nii")
  #imgV2=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V2.nii")
  #imgV3=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_V3.nii")
  
  #V1 = imgV1@.Data
  #V2 = imgV2@.Data
  
  #V1_1=array(V1, c(prod(dim(V1)[1],dim(V1)[2],dim(V1)[3]),3))
  #V2_1=array(V2, c(prod(dim(V2)[1],dim(V2)[2],dim(V2)[3]),3))
  
  #imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/Meditation/con21/data_L1.nii")
  #imgL1=readNIfTI("~/Dropbox/ClusteringDTIonstiefel/RealData/PNC/600085654611/data_L1.nii.gz")
  imgL1=readNIfTI("~/Desktop/new_gen_template.nii.gz")
  
  L1 = imgL1@.Data
  L1_1=array(L1, c(prod(dim(L1)[1],dim(L1)[2],dim(L1)[3]),12))
  L1_1 = L1_1*0
  
  #load('~/Desktop/results_real_6_7_8_1000_beta/outputs/output_nclust_6_real_data_id_12345.RData')
  #load(data_from_result)
  selected_indx = find_validVoxel(cutoff)
  #mcmc_len = length(data_with_init_with_MCMC_samples$MCMC_sample)
  #print(mcmc_len)
  #mcmc_len=210
  #L1_1[selected_indx] = data_with_init_with_MCMC_samples$MCMC_sample[[mcmc_len]]$curr_param$id_arr
  C = dim(soft_cluster_for_data)[2]
  L1_1[selected_indx,] = soft_cluster_for_data
  print(paste0("length of selected data = ",length(L1_1)))
    
  L1_2=array(L1_1, c(dim(L1)[1],dim(L1)[2],dim(L1)[3],12))
  browser()
  
  imgL1@.Data = L1_2
  
  fileName=sprintf('~/Desktop/new_gen_L1_%s',cutoff)
  writeNIfTI(imgL1,fileName)
  
  
}















