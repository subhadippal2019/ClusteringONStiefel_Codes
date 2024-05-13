source("utility.R")
load_src_libs()
source("analyse_MCMC_sample.R")

#### for true cluster = 3
for(i in c(2,3,4,5)){
  for(j in 1:50){  
    fName = sprintf('~/Desktop/DTI_project/simulated_data/results_nclust_3/outputs/output_nclust_%d_data_id_%d_true_nClust_3.RData',i,j)
    load(fName)
    dic = calculate_DIC_4_5_6_7_8(data_with_init_with_MCMC_samples)
    #dic5 = calculate_DIC_5(data_with_init_with_MCMC_samples)
    #dic6 = calculate_DIC_6(data_with_init_with_MCMC_samples)
    
    str1 = sprintf('output_nclust_%d_data_id_%d_true_nClust_3\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',i,j,dic$val4,dic$pd4,dic$val5,dic$pd5,dic$val6,dic$pd6,dic$val7,dic$pd7,dic$val8,dic$pd8)
    
    write.table(str1,file = "./dic_4_5_6_7_8.txt",append=T,row.names=F,col.names=F,quote=F,sep='')
    #cat(str1,"\n")
  }
}

A=read.table("./true_4_dic_4_5_6_7_8.txt")
guess_2_clust_dic_4 = mean(A[1:50,2])
guess_2_clust_dic_5 = mean(A[1:50,4])
guess_2_clust_dic_6 = mean(A[1:50,6])
guess_2_clust_dic_7 = mean(A[1:50,8])
guess_2_clust_dic_8 = mean(A[1:50,10])

guess_3_clust_dic_4 = mean(A[51:100,2])
guess_3_clust_dic_5 = mean(A[51:100,4])
guess_3_clust_dic_6 = mean(A[51:100,6])
guess_3_clust_dic_7 = mean(A[51:100,8])
guess_3_clust_dic_8 = mean(A[51:100,10])

guess_4_clust_dic_4 = mean(A[101:150,2])
guess_4_clust_dic_5 = mean(A[101:150,4])
guess_4_clust_dic_6 = mean(A[101:150,6])
guess_4_clust_dic_7 = mean(A[101:150,8])
guess_4_clust_dic_8 = mean(A[101:150,10])

guess_5_clust_dic_4 = mean(A[151:200,2])
guess_5_clust_dic_5 = mean(A[151:200,4])
guess_5_clust_dic_6 = mean(A[151:200,6])
guess_5_clust_dic_7 = mean(A[151:200,8])
guess_5_clust_dic_8 = mean(A[151:200,10])

guess_6_clust_dic_4 = mean(A[201:250,2])
guess_6_clust_dic_5 = mean(A[201:250,4])
guess_6_clust_dic_6 = mean(A[201:250,6])
guess_6_clust_dic_7 = mean(A[201:250,8])
guess_6_clust_dic_8 = mean(A[201:250,10])


#### for true cluster = 4
for(i in c(6)){
  for(j in 20:50){  
    fName = sprintf('~/Desktop/DTI_project/simulated_data/results_nclust_4/outputs/output_nclust_%d_data_id_%d_true_nClust_4.RData',i,j)
    load(fName)
    dic = calculate_DIC_4_5_6_7_8(data_with_init_with_MCMC_samples)
    #dic5 = calculate_DIC_5(data_with_init_with_MCMC_samples)
    #dic6 = calculate_DIC_6(data_with_init_with_MCMC_samples)
    
    str1 = sprintf('output_nclust_%d_data_id_%d_true_nClust_4\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',i,j,dic$val4,dic$pd4,dic$val5,dic$pd5,dic$val6,dic$pd6,dic$val7,dic$pd7,dic$val8,dic$pd8)
    
    write.table(str1,file = "./true_4_dic4_5_6_7_8.txt",append=T,row.names=F,col.names=F,quote=F,sep='')
  }
}

x = matrix(rep(0,20),ncol=5)
for(i in c(2,3,4,5)){
  for(j in c(4,5,6,7,8)){
    var_name = paste0("guess_",i,"_clust_dic_",j)
    x[i-1,j-3] = get(var_name)
  }
}

for(l in 1:5){
  y = matrix(rep(0,250),ncol=50)
  for(i in c(2,3,4,5,6)){
    k = ((i-2)*50+1):((i-1)*50)
    y[i-1,] = A[k,2*l]
  }

  for(i in c(2,3,4,5,6)){
    print(sum(apply(y,2,which.min)==(i-1)))
    
  }
  print("####")
}

result_NEO = NULL

for (i in 4:8){
  fileName = sprintf("../RealData/NEO/output/output_nclust_%d_realNEO_26094.RData",i)
  tmp1 = strsplit(fileName,"/")[[1]][5]
  str1 = sprintf("%s_%s_%s",strsplit(tmp1,"_")[[1]][2],strsplit(tmp1,"_")[[1]][3],strsplit(tmp1,"_")[[1]][4])
  result_NEO$sample = str1
  load(fileName)
  res = calculate_DIC_5(data_with_MCMC_sample)
  result_NEO$nclust = i
  result_NEO$val = res$val
  result_NEO$pd = res$pd
  write.table(result_NEO,sep="\t",quote=F,row.name=F,col.names=F,append=T,file="DIC_NEO.txt")
}

