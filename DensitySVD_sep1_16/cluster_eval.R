#### cluster evaluation
#### https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
#### April 24, 2017

library(clues)
library(NMI)


##### example:
## load("~/Desktop/DTI_project/simulated_data/results_nclust_3/outputs/output_nclust_3_data_id_17_true_nClust_3.RData")
## load("~/Desktop/DTI_project/simulated_data/ML_dataset_3_17.RData")
## cluster_eval(data_with_init_with_MCMC_samples$MCMC_sample[[2000]]$curr_param$id_arr,L$clust)

f.calculate.tp.tn.fp.fn <- function(ground.truth.membership, clue.membership) {
  if (length(ground.truth.membership) != length(clue.membership)) {
    cat("Vector sizes differ.")
    return()
  }
  
  gt.M <- ground.truth.membership
  cl.M <- clue.membership
  
  contingency <- matrix(0, 2, 2)
  row.names(contingency) <- c("Same Class", "Different Class")
  colnames(contingency) <- c("Same Cluster", "Different Clusters")
  
  result <- data.frame(tp=0, tn=0, fp=0, fn=0)
  
  idxcomb <- combn(1:length(gt.M), 2)
  for (i in 1:ncol(idxcomb)) {
    id1 <- idxcomb[1,i]
    id2 <- idxcomb[2,i]
    if (gt.M[id1] == gt.M[id2] && cl.M[id1] == cl.M[id2]) {
      # TP
      contingency[1,1] = contingency[1,1] + 1
      result$tp <- contingency[1,1]
    } else if (gt.M[id1] != gt.M[id2] && cl.M[id1] != cl.M[id2]) {
      # TN
      contingency[2,2] = contingency[2,2] + 1
      result$tn <- contingency[2,2]
    } else if (gt.M[id1] == gt.M[id2] && cl.M[id1] != cl.M[id2]) {
      # FN
      contingency[1,2] = contingency[1,2] + 1
      result$fn <- contingency[1,2]
    } else if (gt.M[id1] != gt.M[id2] && cl.M[id1] == cl.M[id2]) {
      # FP
      contingency[2,1] = contingency[2,1] + 1
      result$fp <- contingency[2,1]
    }
  }
  
  print(contingency)
  return(result)
}

f.fscore <- function(contingency.values, beta=1) {
  val <- contingency.values
  precision <- val$tp / (val$tp + val$fp)
  recall <- val$tp / (val$tp + val$fn)
  ((beta^2 + 1) * precision * recall) / (beta^2 * precision + recall)
}



purity_calc <- function(cluster, truth) {
  k <- max(cluster, truth)
  cluster <- factor(cluster, levels = 1:k)
  truth <- factor(truth, levels = 1:k)
  m <- length(cluster)
  mi <- table(cluster)
  
  cnts <- split(truth, cluster)
  cnts <- sapply(cnts, FUN = function(n) table(n))
  p <- sweep(cnts, 1, rowSums(cnts), "/")
  p[is.nan(p)] <- 0
  
  res = sum(apply(p, 1, max) * mi/m)
  return(res)
}


cluster_eval <-function(c1,c2){
  
  #c1=c(1,1,1,1,1,2,2,2,2,2,1,3,3,3,3,1,1)
  #c2=c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3)
  
  PUR = purity_calc(c1,c2)  ### purity
  RI = adjustedRand(c1,c2,"Rand") ### Rand Index
  ARI = adjustedRand(c1,c2,"HA")  ### Hubert and Arabie’s adjusted Rand index
  ARI2 = adjustedRand(c1,c2,"MA") ### Morey and Agresti’s adjusted Rand index
  JI = adjustedRand(c1,c2,"Jaccard") ### Jaccard Index
  
  L = length(c1)
  data1 = cbind(seq(1,L),c1)
  data2 = cbind(seq(1,L),c2)
  
  NMI = NMI(data1,data2) ### Normalized mutual information
  
  contingency.values <- f.calculate.tp.tn.fp.fn(c1, c2)
  
  F1 = f.fscore(contingency.values,1)
  F2 = f.fscore(contingency.values, 2)
  F5 = f.fscore(contingency.values, 5)
  F05 = f.fscore(contingency.values, 0.5)
  res = NULL
  
  res$PUR = PUR
  res$RI = unname(RI)
  res$ARI = unname(ARI)
  res$ARI2 = unname(ARI2)
  res$JI = unname(JI)
  res$NMI = NMI$value
  res$F05 = F05
  res$F1 = F1
  res$F2 = F2
  res$F5 = F5
  
  return(res)
}

dirName3 = '~/Desktop/DTI_project/simulated_data/results_nclust_3/outputs'
dirName4 = '~/Desktop/DTI_project/simulated_data/results_nclust_4/outputs'
dirName1 = '~/Desktop/DTI_project/simulated_data'

fid = file("true_3_predicted_3_cluster_eval.txt","w")

str = sprintf('PUR\tRI\tARI\tARI2\tJI\tNMI\tF05\tF1\tF2\tF5\n')
cat(str,file=fid)

for (i in 1:50){
 fName1 = sprintf('%s/ML_dataset_3_%d.RData',dirName1,i) 
 fName3 = sprintf('%s/output_nclust_3_data_id_%d_true_nClust_3.Rdata',dirName3,i)
 load(fName1)
 c1 = L$clust
 load(fName3)
 c2 = data_with_init_with_MCMC_samples$MCMC_sample[[2000]]$curr_param$id_arr
 
 res = cluster_eval(c1,c2)
 
 str = sprintf('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n',res$PUR, res$RI,res$ARI, res$ARI2, res$JI, 
                                                          res$NMI,res$F05, res$F1, res$F2, res$F5)
 
 cat(str,file=fid)
 
 rm(L)
 rm(data_with_init_with_MCMC_samples)
}

fclose(fid)

#####################

fid = file("true_4_predicted_4_cluster_eval.txt","w")
str = sprintf('PUR\tRI\tARI\tARI2\tJI\tNMI\tF05\tF1\tF2\tF5\n')
cat(str,file=fid)

for (i in 1:50){
 fName2 = sprintf('%s/ML_dataset_4_%d.RData',dirName1,i) 
 fName4 = sprintf('%s/output_nclust_4_data_id_%d_true_nClust_4.Rdata',dirName4,i)
 load(fName2)
 c1 = L$clust
 load(fName4)
 c2 = data_with_init_with_MCMC_samples$MCMC_sample[[2000]]$curr_param$id_arr
 
 res = cluster_eval(c1,c2)
 
 str = sprintf('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n',res$PUR, res$RI,res$ARI, res$ARI2, res$JI, 
                                                          res$NMI,res$F05, res$F1, res$F2, res$F5)
 
 cat(str,file=fid)
 
 
 rm(L)
 rm(data_with_init_with_MCMC_samples)
 
}

fclose(fid)


A=read.table("true_3_predicted_3_cluster_eval_NMI.txt",header=T)
m3 = rbind(apply(A,2,summary),sd=apply(A,2,sd))

B=read.table("true_4_predicted_4_cluster_eval_NMI.txt",header=T)
m4 = rbind(apply(B,2,summary),sd=apply(B,2,sd))

write.csv(m3,file="cluster_eval_3_NMI_full_summary.csv")
write.csv(m4,file="cluster_eval_4_NMI_full_summary.csv")

write.csv(round(m3,4),file="cluster_eval_3_NMI_round_summary.csv")
write.csv(round(m4,4),file="cluster_eval_4_NMI_round_summary.csv")
