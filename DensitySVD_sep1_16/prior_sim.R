###### prior simulation
###### 21st May, 2017
###### For 2x2 D matrix in independent prior structure

source("utility.R")


# prior_d_density <- function(y,D,alpha,eta,dimIndex,nu_0,approx=0){
#   
#   ##out=d_density_d1_log(y=(1:4000)/2000,D=c(0.7,5.4),S=c(14,15),1,10,1)
#   ##plot(exp(out-max(out)))
# 
#   val=0*seq(1:length(y));
#   
#   for (i in  1:length(y)){
#     
#     x=y[i]
#     if(x > 100){
#       approx = 1
#     }
#     D[dimIndex]=x
#     eigenValues=D^2/4
#     d = 0.0
#     
#     if(approx==1){
#       log_hyper0F1_val = approx_hyper_log(eigenValues)  #### defined in utility.R 
#       #print(paste0("******######******** ",hyper0F1_val))
#     }else{
#       hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
#       log_hyper0F1_val = log(hyper0F1_val)
#     }
#     val[i]=nu_0*(x*eta[dimIndex])- nu_0 * log_hyper0F1_val
#       
#   }
#   #write.table(val,file="./density1.txt")
#   return(val)  
# }


# prior_d_bivariate_density <- function(y1,y2,D,alpha,eta,dimIndex,nu_0,approx=0){
#   
#   ##out=d_density_d1_log(y=(1:4000)/2000,D=c(0.7,5.4),S=c(14,15),1,10,1)
#   ##plot(exp(out-max(out)))
#   
#   r1 = length(y1)
#   c1 = length(y2)
#   val = matrix(rep(0,r1*c1),nrow=r1)
#   
#   for(i in 1:length(y1)){
#     for (j in  1:length(y2)){
#     
#       x = c(y1[i],y2[j])
#       
#       if(any(x > 100)){
#         approx = 1
#       }
#       D = x
#       eigenValues=D^2/4
#       d = 0.0
#     
#       if(approx==1){
#         log_hyper0F1_val = approx_hyper_log(eigenValues)  #### defined in utility.R 
#         #print(paste0("******######******** ",hyper0F1_val))
#       }else{
#         hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,d)[[3]]
#         log_hyper0F1_val = log(hyper0F1_val)
#       }
#     
#       val[i,j]=nu_0*sum(x*eta)- nu_0 * log_hyper0F1_val
#     }
#   }
#   #write.table(val,file="./density1.txt")
#   return(val)  
# }


prior_d1_d2_bivariate_density <- function(y1,y2,eta,nu_0,approx=0){
  
  r1 = length(y1)
  c1 = length(y2)
  val = matrix(rep(0,r1*c1),nrow=r1)
  
  for(i in 1:length(y1)){
    for (j in  1:length(y2)){
      
      x = c(y1[i],y2[j])
      
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
  #write.table(val,file="./density1.txt")
  return(val)  
}


load_src_libs()
#### parameters of our prior
#alpha = c(1,1)
nu_0 = 10
eta = c(0.9,0.8)#c(0.9166642,0.8802722)

#D = NULL
#dimIndex = 1


### univariate
#y = seq(0.1,10,by=0.01)
#out = prior_d_density(y1,y2,D,alpha,eta,dimIndex,nu_0,approx=0)
#plot(y,exp(out-max(out)),type='l')

fd <- file("./density1.txt", "w")

### bivariate
y1 = seq(1,20,by=0.5)
y2 = seq(1,10,by=0.5)
out = prior_d1_d2_bivariate_density(y1,y2,eta,nu_0,approx=0)
#pdf("bivariate_conjugate_prior.pdf")
exp_out = exp(out-max(out))
#image(exp_out)
#dev.off()

for(i in 1:length(y1)){
  for (j in  1:length(y2)){
    cat(y1[i],"\t",y2[j],"\t",exp_out[i,j],"\n",file=fd,append=TRUE)
  }
}

close(fd)

library(ggplot2)
library(RColorBrewer)
library(ggthemes)

data = read.table('density1.txt')

x=data[,1]
y=data[,2]
z=data[,3]
# build your data.frame
df <- data.frame(Eigenvalue_1=x, Eigenvalue_2=y, Prior_density=z)

# build color Palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

ggplot(df) + 
  aes(x = Eigenvalue_1, y = Eigenvalue_2, z = Prior_density, fill = Prior_density) + 
  geom_tile() + 
  ggtitle("Prior density plot") +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_distiller(palette="Spectral", na.value="white") + 
  theme(legend.position="bottom")

#require(plot3D)
#persp3D(z = exp_out, theta = 30)

#image2D(exp_out)
#par(new=TRUE)
#contour2D(exp_out,xlim=c(0,1),ylim=c(0,1))


