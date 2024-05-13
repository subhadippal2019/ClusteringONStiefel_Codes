#### rgl to create density / contours on 3D sphere

library(rgl)
library(reshape)

data(volcano) 
y <- 2 * volcano # Exaggerate the relief
x <- 10 * (1:nrow(y)) # 10 meter spacing (S to N)
z <- 10 * (1:ncol(y)) # 10 meter spacing (E to W)
ylim <- range(y)
ylen <- ylim[2] - ylim[1] + 1
colorlut <- terrain.colors(ylen) # height color lookup table
col <- colorlut[ y-ylim[1]+1 ] # assign colors to heights
rgl.open()
rgl.surface(x, z, y, color=col, back="lines")

rgl.clear(type = c("shapes"))
calls = read.delim('call_data.tsv', header = T)

bin_size = 0.18
calls$long_bin = cut(calls$long, seq(min(calls$long), max(calls$long), bin_size))
calls$lat_bin = cut(calls$lat, seq(min(calls$lat), max(calls$lat), bin_size))
calls$total = log(calls$total) / 3 #flatten out totals

calls = melt(calls[,3:5])
calls = cast(calls, lat_bin~long_bin, fun = sum, fill = 0)
calls = calls[,2:(ncol(calls)-1)]
calls = as.matrix(calls)

x = (1: nrow(calls))
z = (1: ncol(calls))
rgl.surface(x, z, calls)
rgl.bringtotop()

rgl.pop()
# nicer colored plot
ylim <- range(calls)
ylen <- ylim[2] - ylim[1] + 1
col <- topo.colors(ylen)[ calls-ylim[1]+1 ]
x = (1: nrow(calls))
z = (1: ncol(calls))

rgl.bg(sphere=FALSE, color=c("black"), lit=FALSE)
rgl.viewpoint( theta = 300, phi = 30, fov = 170, zoom = 0.03)
rgl.surface(x, z, calls, color = col, shininess = 10)
rgl.bringtotop()

# rgl demo: rgl-bivar.r
# author: Daniel Adler
# $Id$

rgl.demo.bivar <- function()
{
  require(MASS);
  
  # parameters:
  n<-50; ngrid<-40
  
  # generate samples:
  set.seed(31415)
  x<-rnorm(n); y<-rnorm(n)
  
  # estimate non-parameteric density surface via kernel smoothing
  denobj<-kde2d(x, y, n=ngrid)
  den.z <-denobj$z
  
  # generate parametric density surface of a bivariate normal distribution
  xgrid <- denobj$x
  ygrid <- denobj$y
  bi.z <- dnorm(xgrid)%*%t(dnorm(ygrid))
  
  # visualize:
  zscale<-20
  
  # New window
  open3d()
  
  # clear scene:
  clear3d("all")
  
  # setup env:
  bg3d(color="#887777")
  light3d()
  
  # Draws the simulated data as spheres on the baseline
  spheres3d(x,y,rep(0,n),radius=0.1,color="#CCCCFF")
  
  # Draws non-parametric density
  surface3d(xgrid,ygrid,den.z*zscale,color="#FF2222",alpha=0.5)
  
  # Draws parametric density
  surface3d(xgrid,ygrid,bi.z*zscale,color="#CCCCFF",front="lines") 
}

rgl.demo.bivar()
