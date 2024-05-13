data <- read.csv('~/Dropbox/ClusteringDTIonstiefel/RealData/NEO/NEO1.csv', sep=',', stringsAsFactors = T,header=T)
L.col = 5
t.col = 6
nu.col = 7

L =  data[,L.col]*(pi/180)
t =  data[,t.col]*(pi/180)
nu = data[,nu.col]*(pi/180)

r = sqrt(sin(t)*sin(t) + cos(t)*cos(t)*sin(nu-L)*sin(nu-L))
x_1 = cbind(cos(t)*cos(L),cos(t)*sin(L),sin(t))
x_2 = cbind(sin(t)*sin(nu),-sin(t)*cos(nu),-cos(t)*sin(nu-L))/r

N = length(r)
x = array(rep(0,3*2*N), c(3, 2, N))
for (i in 1:N){
  x[,,i] = cbind(t(x_1[i,]),t(x_2[i,]))
}

neo_192 = x
save(neo_192,file="neo_192_data.RData")