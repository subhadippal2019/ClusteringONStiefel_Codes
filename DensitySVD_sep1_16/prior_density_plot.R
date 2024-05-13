###### prior density plot function
###### 7th July, 2017
###### For 2x2 D matrix in independent prior structure

library(ggplot2)
library(RColorBrewer)
library(ggthemes)

prior_d1_d2_bivariate_density <- function(d1,d2,eta,nu_0,approx=0){
  
  r1 = length(d1)
  c1 = length(d2)
  val = matrix(rep(0,r1*c1),nrow=r1)
  
  for(i in 1:length(d1)){
    for (j in  1:length(d2)){
      
      x = c(d1[i],d2[j])
      
      if(any(x > 100)){
        approx = 1
      }
      
      D = x
      eigenValues=D^2/4
      d = 0.0
      
      if(approx==1){
        log_hyper0F1_val = approx_hyper_log(eigenValues)  #### defined in utility.R 
        #print(paste0("******######******** ",hyper0F1_val))
      }else{
        hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
        log_hyper0F1_val = log(hyper0F1_val)
      }
      
      val[i,j]=nu_0*sum(x*eta) - nu_0 * log_hyper0F1_val
      
    }
  }
  return(val)  
}


plot_prior_density <- function(nu_0,eta){

  source("utility.R")
  load_src_libs()
  
  d1min = 0.05
  d1max = 18
  d1by = 0.05
  
  d2min = 0.05 
  d2max = 18
  d2by = 0.05
  
  ### bivariate distribution
  d1 = seq(d1min,d1max,by=d1by)
  d2 = seq(d2min,d2max,by=d2by)
  
  out = prior_d1_d2_bivariate_density(d1,d2,eta,nu_0,approx=0)
  
  ## normalize such that max is 1.0
  exp_out = exp(out-max(out))
  
  #fd <- file("./density1.txt", "w")
  
  L1 = length(d1)
  L2 = length(d2)
  data = matrix(rep(0,L1*L2*3),ncol=3)
  
  k = 0
  for(i in 1:length(d1)){
    for (j in  1:length(d2)){
      k = k+1
      data[k,] = c(d1[i],d2[j],exp_out[i,j])
      #cat(y1[i],"\t",y2[j],"\t",exp_out[i,j],"\n",file=fd,append=TRUE)
    }
  }
  
  #close(fd)
  
  #data = read.table('prior_density.txt')
  
  x = data[,1]
  y = data[,2]
  z = data[,3]
  
  
  # build the data.frame
  df = data.frame(d_1=x, d_2=y, Prior_density=z)
  
  
  mode = which(exp_out==max(exp_out),arr.ind = T)
  mode_x = d1[mode[1]]
  mode_y = d2[mode[2]]
  

  ggplot(df) + 
    aes(x = d_1, y = d_2, z = Prior_density, fill = Prior_density) + 
    geom_tile() +
    geom_contour(color = "black", alpha = 0.5) +
    scale_fill_distiller(palette="Spectral", na.value="white") + 
    theme(legend.position="none",text = element_text(size=11)) +
    geom_segment(data=df,mapping=aes(x=mode_x, y=0, xend=mode_x, yend=mode_y), linetype = 2, size=0.4, color="black") + 
    geom_segment(data=df,mapping=aes(x=0, y=mode_y, xend=mode_x, yend=mode_y), linetype = 2, size=0.4, color="black") +
    scale_x_continuous( breaks = c(0,mode_x,10,15,18), labels = c(0,mode_x,10,15,18)) + 
    scale_y_continuous( breaks = c(0,mode_y,10,15,18), labels = c(0,mode_y,10,15,18)) + 
    labs(x=expression(d[1])) +   
    labs(y=expression(d[2])) +
    
    theme(axis.title.x = element_text(size = rel(1.1), angle = 00)) + 
    theme(axis.title.y = element_text(size = rel(1.1), angle = 90)) 

} 


##### from given D_mode find out the appropriate \eta
find_eta_for_mode <-function(d1_m,d2_m){
  
  D = c(d1_m,d2_m)
  eigenValues=D^2/4
  d = 0.0
  hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
  partial_d1_hyper0F1_val = .C("partial_d1_hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
  partial_d2_hyper0F1_val = .C("partial_d2_hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
  
  eta_m = NULL
  eta_m[1] = partial_d1_hyper0F1_val / hyper0F1_val
  eta_m[2] = partial_d2_hyper0F1_val / hyper0F1_val
  
  return(eta_m)
}


### main function
#main_plot <-function(){

  source("utility.R")
  load_src_libs()
  #### parameters of our prior
  #nu_0_arr = c(8,14,20,30)
  eta_mat = matrix(c(0.89,0.5, 0.5,0.89, 0.94,0.94),ncol=2,byrow=T) 
  ##############
  
  ### nu_0 = 8,14,20
  ### eta = [0.89,0.5], [0.5,0.89], [0.8,0.8], [0.94,0.94] with 35
  
  #for(nu_0 in nu_0_arr){
    #for(j in c(1,2,3) ){
    #for(j in c(1) ){
        
      nu_0 = 15
  
      d1_m = 7
      d2_m = 5
  
      #eta = eta_mat[2,]
      eta = eta_m = find_eta_for_mode(d1_m,d2_m)
      
      #pngname = sprintf("/Users/subhajit/Dropbox/ClusteringDTIonstiefel/DensitySVD_sep1_16/eps_plots_1/plot_nu0_%d_eta_%0.2f_%0.2f.png",nu_0,eta[1],eta[2])
      #cat("png ", j, ": ",pngname," will be generated\n")
      pngname = sprintf("/Users/subhajit/Dropbox/ClusteringDTIonstiefel/DensitySVD_sep1_16/eps_plots_1/plot_nu0_%d_mode_7_5.png",nu_0,eta[1],eta[2])
      png(pngname,res=200,width=800,height=800)
      #pdf(pngname)
      plot_prior_density(nu_0,eta)
      #dev.off()
      cat("eta_m is ", eta, "\n")
      
    #}
    
  #}
      
      

#}