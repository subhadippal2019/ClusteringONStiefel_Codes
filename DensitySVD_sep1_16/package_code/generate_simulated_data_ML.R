generate_simulated_data_ML<-function(nClust,run_id,N=100){

  if((nClust != 3) && (nClust !=4)){
    stop("only 3 or 4 clusters!!")
  }

  #source("bivariate_NR_method_d1_d2.R")
  #source("stiefel_SVD.R")
  ###################
  ### Description:
  
  ### Simulated data generation for fixed number of clusters 
  ### Visualize in some random direction on 2D
  ### Clustering heatmap visualization
  
  ###################
  
  ##################
  #library(rstiefel)
  ##################
  
  if(nClust == 3){
    M_hyper_arr = array(c(1,0,0,0,0,1,0,1,0,0,0,-1,0,0,1,1,0,0), c(3, 2, nClust))
  }
  if(nClust == 4){
    M_hyper_arr = array(c(1,0,0,0,0,-1,0,-1,0,0,0,1,0,0,1,0,-1,0,0,1,0,-1,0,0), c(3, 2, nClust))
  }
  
  #M_hyper_arr = array(c(1,0,0,0,1,0,1,0,0,0,-1,0,1/sqrt(2),0,1/sqrt(2),0,1,0), c(3, 2, 3))  
  
  ### mean matrix generation nClust matrices
  M = array(rep(0,3*2*nClust), c(3, 2, nClust))
  for(i in 1:nClust){
    #M[,,i] = rmf.matrix(M_hyper_arr[,,i]) 
    M[,,i] = M_hyper_arr[,,i]
  }
  
  D = array(rep(0,2*2*nClust), c(2, 2, nClust))
  for(i in 1:nClust){
    s = 0.5*sort(rgamma(c(1,1),c(160,100),c(20,20)),decreasing=TRUE)
    D[,,i] = diag(s+i)
  }
 
  V = array(rep(0,2*2*nClust), c(2, 2, nClust))
  for(i in 1:nClust){
    V[,,i] = matrix(c(1,0,0,1))
  }
  ###### main parameter
  F = array(rep(0,3*2*nClust), c(3, 2, nClust))
  
  ##########################
  ### to keep the sign thing consistent
  
  for(i in 1:nClust){
    M_V_list = change_sign_M_V(M[,,i],V[,,i])
    M[,,i] = M_V_list$M
    V[,,i] = M_V_list$V
    F[,,i] = M[,,i]%*%D[,,i]%*%V[,,i]
  }
  
  
  ##########################
  data = array(rep(0,3*2*N), c(3, 2, N)); 
  ## 0.2 0.4 0.4
  cluster_for_data = rep(0,N)
  prob_vec = rep(1,nClust)/nClust
  for(i in 1:N){
    r = rmultinom(1,1,prob_vec)
    id = which(r==1)
    F1 = M[,,id]%*%D[,,id]%*%V[,,id]
    
    data[,,i] = rmf.matrix(F1) 
    cluster_for_data[i] = id
  }
  
  #print(data)
  ### reindexing the data to store the data according to the cluster order
  B = sort(cluster_for_data,index.return	= TRUE)
  L = NULL
  L = list(data=data[,,B$ix],clust=cluster_for_data[B$ix],curr_param=list(M=M,D=D,V=V,F=F))
  #print(L)
  
  file_name = sprintf("ML_dataset_%d_%d.RData",nClust,run_id)
  save(L,file=file_name)
  #return(L) 
}

###### some random projection of data on 2D
visualize_data<-function(L){
  
  data = L$data
  r1 = runif(dim(data)[1],-1,1)
  random_direction_vec = r1/sqrt(sum(r1*r1))
  
  N = dim(data)[3]
  proj = matrix(rep(0,2*N), c(2, N)); 
  for(i in 1:N){
      proj[,i] = random_direction_vec%*%data[,,i]   
  }
  arr = c('red','blue','green')
  plot(proj[1,],proj[2,],col=arr[L$clust])
  
}


####### clustering visualization
calc_dist_data<-function(L){
  data = L$data
  N = dim(data)[3]
  
  data1 = array(rep(0,3*2*N), c(3, 2, N));
  nClust = length(unique(L$clust))
  
#   ID1 = NULL
#   ID2 = NULL
#   ID3 = NULL
#   cnt1 = 0
#   cnt2 = 0
#   cnt3 = 0
#   
#   
#   
#   for(i in 1:N){
#     if(L$clust[i] == 1){
#       cnt1 = cnt1+1
#       ID1[cnt1] = i
#     }
#     if(L$clust[i] == 2){
#       cnt2 = cnt2+1
#       ID2[cnt2] = i
#     }
#     if(L$clust[i] == 3){
#       cnt3 = cnt3+1
#       ID3[cnt3] = i
#     }
#     
#   }
  
  id_start = 1
  for(i in 1:nClust){
    cnt1 = sum(L$clust == i)
    id_end =  id_start + cnt1-1
    data1[,,id_start:id_end] = data[,,which(L$clust == i)]
    id_start = id_end+1
  }
  
  #data1[,,(cnt1+1):(cnt1+cnt2)] = data[,,ID2]
  #data1[,,(cnt1+cnt2+1):N] = data[,,ID3]
  
  
  dist_mat = matrix(rep(0.0,N*N),nrow=N,ncol=N)
  
  for(i in 1:N){
    for(j in 1:N){
      
      d = 2.0-sum(diag(t(data1[,,i])%*%data1[,,j]))
      dist_mat[i,j] = d
    }
  }
  
  #save(dist_mat,file="dist_mat.Rdata")
  heatmap(dist_mat,Rowv=NA,Colv=NA)
}


  
