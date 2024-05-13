#### co-clustering matrix
# x, y, z : numeric vectors corresponding to
#  the coordinates of points
# axis.col : axis colors
# xlab, ylab, zlab: axis labels
# show.plane : add axis planes
# show.bbox : add the bounding box decoration
# bbox.col: the bounding box colors. The first color is the
# the background color; the second color is the color of tick marks
rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{ 
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)
  
  # Add plane
  if(show.plane) 
    xlim <- xlim/1.1; zlim <- zlim /1.1
  rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
             z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
  
  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 3, ylen = 3, zlen = 3) 
  }
}




rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

load("./NEC_results_july29/output/july27_output_nclust_4_realNEC.RData")

L = length(data_with_MCMC_sample$MCMC_sampl)
N = length(data_with_MCMC_sample$MCMC_sample[[1]]$curr_param$id_arr)

co_cluster_mat = matrix(rep(0,N*N),ncol=N)

for (i in 1:N){
  for(j in i:N){
    for(k in 1:L){
      elem1 = data_with_MCMC_sample$MCMC_sample[[k]]$curr_param$id_arr[i]
      elem2 = data_with_MCMC_sample$MCMC_sample[[k]]$curr_param$id_arr[j]
      if(elem1 == elem2)
        co_cluster_mat[i,j] = co_cluster_mat[i,j] + 1     
    }
  }
}

for (i in 1:N){
  for(j in i:N){
    co_cluster_mat[j,i] = co_cluster_mat[i,j]  
  }
}

for (i in 1:N){
  for(j in 1:N){
    co_cluster_mat[i,j] = co_cluster_mat[i,j] / L 
  }
}

png("co_clustering_prob.png",width=480,height=540)
library(reshape2)
library(ggplot2)
ggplot(melt(co_cluster_mat), aes(Var1,Var2, fill=value)) + geom_raster() + labs(x="Index of NEC object",y="Index of NEC object") + theme(legend.position="bottom") + 
  guides(fill =
           guide_legend(
             title = "Co-occurrence probability  ",
             title.theme = element_text(
               size = 15,
               face = "italic",
               colour = "black",
               angle = 0
             )
           )
  )
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist = N^2
iter_indx = 0

for(k in 1:L){
  tmp_mat = matrix(rep(0,N*N),ncol=N)
  for (i in 1:N){
    for(j in i:N){
      elem1 = data_with_MCMC_sample$MCMC_sample[[k]]$curr_param$id_arr[i]
      elem2 = data_with_MCMC_sample$MCMC_sample[[k]]$curr_param$id_arr[j]
      if(elem1 == elem2)
        tmp_mat[i,j] = 1 
    }
  }
  dist_tmp = sum(abs(tmp_mat-co_cluster_mat))
  if(dist_tmp < dist){
    dist = dist_tmp
    iter_indx = k   ### 593
  }
}

#### optimal cluster assignment
idx = data_with_MCMC_sample$MCMC_sample[[iter_indx]]$curr_param$id_arr

#### plotting points on 3D sphere
#N = 10
x1 = rep(0,N)
y1 = rep(0,N)
z1 = rep(0,N)
t1 = rep(0,N)

x2 = rep(0,N)
y2 = rep(0,N)
z2 = rep(0,N)
t2 = rep(0,N)


for (i in 1:N){
  
  x1[i] = data_with_MCMC_sample$data[1,1,i]
  y1[i] = data_with_MCMC_sample$data[2,1,i]
  z1[i] = data_with_MCMC_sample$data[3,1,i]
  t1[i] = sqrt(x1[i]*x1[i]+y1[i]*y1[i]+z1[i]*z1[i])
  x1[i] = x1[i]/t1[i]
  y1[i] = y1[i]/t1[i]
  z1[i] = z1[i]/t1[i]
  
  
  x2[i] = data_with_MCMC_sample$data[1,2,i]
  y2[i] = data_with_MCMC_sample$data[2,2,i]
  z2[i] = data_with_MCMC_sample$data[3,2,i]
  t2[i] = sqrt(x2[i]*x2[i]+y2[i]*y2[i]+z2[i]*z2[i])
  x2[i] = x2[i]/t2[i]
  y2[i] = y2[i]/t2[i]
  z2[i] = z2[i]/t2[i]
  
}

#plot3d(x1,y1,z1)


library(RColorBrewer)
library(rgl)

CC = 4

bbox.col = c("#333377","black")
axis.col = "black"
open3d()
rgl_init()
abclines3d(0, 0, 0, a = 2*diag(3), col = "black",lwd  = 3.0)
#rgl.spheres(0,0,0,radius=1,color=c("white"),alpha=0.65)
#rgl_add_axes(x1, y1, z1,show.bbox = FALSE,show.bbox = TRUE)
# Compute and draw the ellipse of concentration
sigma =  matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3)
ellips <- ellipse3d(sigma, centre=c(0, 0, 0),level=0.205)
shade3d(ellips, col = "#D95F02", alpha = 0.1, lit = FALSE)
wire3d(ellips, col = "#D95F02",  lit = FALSE)
#aspect3d(1,1,1)

for (c in 1:CC){
  indx_tmp = which(idx == c)
  spheres3d(x1[indx_tmp],y1[indx_tmp],z1[indx_tmp],col=c,radius=0.03)
  #rgl.points(x1[indx_tmp],y1[indx_tmp],z1[indx_tmp],color = c,size = 5.0)
}
#rgl.postscript("eig1.pdf","pdf")
#rgl.close()

rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
         emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
         xlen = 3, ylen = 3, zlen = 3) 

#text3d(matrix(c(0,1.2,0,1.2,0,0,0,0,1.2),ncol=3),texts=c('y', 'x', 'z'),adj=3,font=3,cex=2)
text3d(matrix(c(0,0,1.2,0,1.2,0,1.2,0,0),ncol=3),texts=c('z', 'y', 'x'),adj=3,font=3,cex=2)

#### 2nd Eigenvector
open3d()
rgl_init()
#rgl.spheres(0,0,0,radius=1,color=c("white"),alpha=0.65)
abclines3d(0, 0, 0, a = 2*diag(3), col = "black",lwd  = 3.0)
#rgl_add_axes(x1, y1, z1,show.bbox = TRUE,show.plane = FALSE)
# Compute and draw the ellipse of concentration
sigma =  matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3)
ellips <- ellipse3d(sigma, centre=c(0, 0, 0),level=0.205)
shade3d(ellips, col = "#D95F02", alpha = 0.1, lit = FALSE)
wire3d(ellips, col = "#D95F02",  lit = FALSE)
aspect3d(1,1,1)
#rgl.spheres(0,0,0,radius=1,color=c("white"),alpha=0.65)
for (c in 1:CC){
  indx_tmp = which(idx == c)
  spheres3d(x2[indx_tmp],y2[indx_tmp],z2[indx_tmp],col=c,radius=0.03)
  #rgl.points(x2[indx_tmp],y2[indx_tmp],z2[indx_tmp],color = c,size = 5.0)
}
#rgl_add_axes(x1, y1, z1)
#rgl.postscript("eig2.pdf","pdf")
rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
         emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
         xlen = 3, ylen = 3, zlen = 3) 
text3d(matrix(c(0,0,1.2,0,1.2,0,1.2,0,0),ncol=3),texts=c('z', 'y', 'x'),adj=3,font=3,cex=2)

#### example

# set.seed(101)
# n <- 50
# theta <- runif(n,0,2*pi)
# u <- runif(n,-1,1)
# x <- sqrt(1-u^2)*cos(theta)
# y <- sqrt(1-u^2)*sin(theta)
# z <- u
# spheres3d(x[1:20],y[1:20],z[1:20],col="red",radius=0.02,color = rainbow(10))
# spheres3d(x[30:40],y[30:40],z[30:40],col="red",radius=0.02)
# 



##############
bbox.col = c("#333377","black")
axis.col = "black"
##### simulation results
load("~/Desktop/DTI_project/simulated_data/results_nclust_4/outputs/output_nclust_4_data_id_1_true_nClust_4.RData")

L = length(data_with_init_with_MCMC_samples$MCMC_sampl)
N = length(data_with_init_with_MCMC_samples$MCMC_sample[[1]]$curr_param$id_arr)

iter_indx = L
idx = data_with_init_with_MCMC_samples$MCMC_sample[[iter_indx]]$curr_param$id_arr

x1 = rep(0,N)
y1 = rep(0,N)
z1 = rep(0,N)
t1 = rep(0,N)

x2 = rep(0,N)
y2 = rep(0,N)
z2 = rep(0,N)
t2 = rep(0,N)

for (i in 1:N){
  
  x1[i] = data_with_init_with_MCMC_samples$data[1,1,i]
  y1[i] = data_with_init_with_MCMC_samples$data[2,1,i]
  z1[i] = data_with_init_with_MCMC_samples$data[3,1,i]
  t1[i] = sqrt(x1[i]*x1[i]+y1[i]*y1[i]+z1[i]*z1[i])
  x1[i] = x1[i]/t1[i]
  y1[i] = y1[i]/t1[i]
  z1[i] = z1[i]/t1[i]
  
  
  x2[i] = data_with_init_with_MCMC_samples$data[1,2,i]
  y2[i] = data_with_init_with_MCMC_samples$data[2,2,i]
  z2[i] = data_with_init_with_MCMC_samples$data[3,2,i]
  t2[i] = sqrt(x2[i]*x2[i]+y2[i]*y2[i]+z2[i]*z2[i])
  x2[i] = x2[i]/t2[i]
  y2[i] = y2[i]/t2[i]
  z2[i] = z2[i]/t2[i]
  
}
#### 1st Eigenvector of simulation data
open3d()
rgl_init()
abclines3d(0, 0, 0, a = 2*diag(3), col = "black",lwd  = 3.0)
# Compute and draw the ellipse of concentration
sigma =  matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3)
ellips <- ellipse3d(sigma, centre=c(0, 0, 0),level=0.205)
shade3d(ellips, col = "#D95F02", alpha = 0.1, lit = FALSE)
wire3d(ellips, col = "#D95F02",  lit = FALSE)
aspect3d(1,1,1)

CC = 4

for (c in 1:CC){
  indx_tmp = which(idx == c)
  spheres3d(x1[indx_tmp],y1[indx_tmp],z1[indx_tmp],col=c,radius=0.03)
}
rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
         emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
         xlen = 3, ylen = 3, zlen = 3) 
text3d(matrix(c(0,0,1.2,0,1.2,0,1.2,0,0),ncol=3),texts=c('z', 'y', 'x'),adj=3,font=3,cex=2)


#### 2nd Eigenvector of simulation data
open3d()
rgl_init()
abclines3d(0, 0, 0, a = 2*diag(3), col = "black",lwd  = 3.0)
# Compute and draw the ellipse of concentration
sigma =  matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3)
ellips <- ellipse3d(sigma, centre=c(0, 0, 0),level=0.205)
shade3d(ellips, col = "#D95F02", alpha = 0.1, lit = FALSE)
wire3d(ellips, col = "#D95F02",  lit = FALSE)
aspect3d(1,1,1)
for (c in 1:CC){
  indx_tmp = which(idx == c)
  spheres3d(x2[indx_tmp],y2[indx_tmp],z2[indx_tmp],col=c,radius=0.03)
}
rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
         emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
         xlen = 3, ylen = 3, zlen = 3) 
text3d(matrix(c(0,0,1.2,0,1.2,0,1.2,0,0),ncol=3),texts=c('z', 'y', 'x'),adj=3,font=3,cex=2)



