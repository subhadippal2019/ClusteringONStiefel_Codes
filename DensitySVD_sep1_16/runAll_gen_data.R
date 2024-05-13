
##### generating all data for true_nclust = 3 or 4
##### in each case generate 50 (differ by run_id) datasets

#R CMD BATCH --no-save --no-restore '--args ./data1 4' ./runAll_gen_data.R ./run1.log

rm(list=ls(all=TRUE))
args <- commandArgs(TRUE)
source("utility.R")

load_src_libs()
N_tot = as.numeric(args[2])
nClust = 3

data_dir = args[1]

for(i in 1:N_tot){
	generate_simulated_data_ML(nClust,run_id=i,N=400,data_dir)
}
nClust = 4
for(i in 1:N_tot){
	generate_simulated_data_ML(nClust,run_id=i,N=500,data_dir)
}
