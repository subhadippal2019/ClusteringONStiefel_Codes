load("../RealData/NEC/nec_175_data.RData")

N = dim(neo_175)[3]
p = dim(neo_175)[2]
  
dist_nec_175 = matrix(rep(0,175*175),ncol=175)
for(i in 1:N){
  for(j in 1:N){
    dist_nec_175[i,j] = p - sum(diag(t(neo_175[,,i])%*%neo_175[,,j]))
  }
}
nec_175.dist = as.dist(dist_nec_175)

save(nec_175.dist, file="dist_nec_175.RData")
library("rgl")
data.mds = cmdscale(nec_175.dist, k=3) 

#Create x,y refs 
data.x <- (data.mds[,1])
data.y <- (data.mds[,2])
data.z <- (data.mds[,3])

load("../RealData/NEC/output/output_nclust_3_realNEC.RData")
data.clust = data_with_MCMC_sample$MCMC_sample[[1600]]$curr_param$id_arr  

#Plot 
#plot(data.x, data.y,col=as.integer(data.clust))

plot3d(data.x, data.y, data.z,col=as.integer(data.clust),radius=5,lwd=2)
#play3d(spin3d(axis=c(0,1,1), rpm=3), duration=30)