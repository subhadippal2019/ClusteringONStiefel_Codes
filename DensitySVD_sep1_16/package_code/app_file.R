source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\NR_method_d1.R')
source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\analyse_MCMC_sample.R')
source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\bivariate_NR_method_d1_d2.R')
source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\cluster_param_update.R')
source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\finiteMixtureML.R')
source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\generate_simulated_data_ML.R')

source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\preprocess_data.R')
source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\rdensity_d1_d2.R')
source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\stiefel_SVD.R')
#source('C:\\Users\\subha\\Desktop\\manifoldCluster\\package_code\\utility.R')
library(Rcpp)
library(gsl)
sourceCpp('C:/Users/subha/Dropbox/projects/ClusteringDTIonStiefel/DensitySVD_sep1_16/package_code/src/functions_hyper_2by2_CPP.cpp')


load("ML_dataset_3_1.RData")
Data=L$data



#DTI DATA
load("C:\\Users\\subha\\Dropbox\\projects\\ClusteringDTIonStiefel\\RealData\\PNC\\600085654611\\real_DTI_data_0.31_PNC_large.RData")
Data1=real_DTI_data[, , 1:4000]
mcmc=finiteMixtureML(data = Data1, nc = 10, max_iter = 10, run_id = 50)
