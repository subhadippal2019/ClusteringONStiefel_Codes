###########
### utility functions

#' It is a debug print function
#' @param str = input string 
#' @param debug_flag = boolean flag to print or not
#' @return NULL  
print_debug <-function(str,debug_flag){
  if(debug_flag == 1){
    print(str)
  }
}

err_msg <-function(str){
  print(paste0('***********ERROR*************'))
  print(str)
  print(paste0('***********ERROR*************'))
}  


load_src_libs <-function()
{
  require(rstiefel)
  require(gtools)
  require(expm)
  require(MASS)

  source("cluster_param_update.R")
  source("rdensity_d1_d2.R")
  source("NR_method_d1.R")
  source("analyse_MCMC_sample.R")
  source("bivariate_NR_method_d1_d2.R")
  source("stiefel_SVD.R")
  #source("utility.R")
  source("preprocess_data.R")
  source("generate_simulated_data_ML.R")
  source("finiteMixtureML.R")

  dyn.load("./src/functions_hyper_2by2_R.so")	  
}


#### DIC calculation function
