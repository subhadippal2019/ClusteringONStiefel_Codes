initial_cluster_center <- function(){
  library(flexclust)
  load("ML_dataset.Rdata")
  vector_data = t(apply(L$data,3,function(x) {return(as.vector(x))}))
  
  nc = 3
  kf = kccaFamily(which = "angle")
  cl <- kcca(x=vector_data[1:100,], k=nc , family=kf)

  #dist=dist_func
}

dist_func <- function (x, centers) 
{
   z <- matrix(0, nrow(x), ncol = nrow(centers))
   for (k in 1:nrow(centers)) {
     z[, k] <- 2 - sum(x*centers[k, ])
   }
  
}


