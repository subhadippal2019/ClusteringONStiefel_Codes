
#nohup R CMD BATCH --no-save --no-restore '--args ./new_data/ML_dataset_3_1.Rdata 2' ./runAllCode.R ./run_2_1_3.log&

#output_nclust_3_data_id_1_true_nClust_3.RData
# other two files (MCMC and init) have the same naming convention

fid = open("input.txt",'w')

data_dir = "./data1"
Nid = 2
nclustMax = 5
true_nclust = 3
for run_id in range(1,Nid):
	for nclust in range(2,nclustMax):
		str = 'nohup R CMD BATCH --no-save --no-restore \'--args %s/ML_dataset_%d_%d.Rdata %d\' ./runAllCode.R ./run_%d_%d_%d.log&'%(data_dir,true_nclust,run_id,nclust,nclust,run_id,true_nclust)	
		fid.write(str+'\n')

true_nclust = 4
for run_id in range(1,Nid):
	for nclust in range(2,nclustMax):
		str = 'nohup R CMD BATCH --no-save --no-restore \'--args %s/ML_dataset_%d_%d.Rdata %d\' ./runAllCode.R ./run_%d_%d_%d.log&'%(data_dir,true_nclust,run_id,nclust,nclust,run_id,true_nclust)	
		fid.write(str+'\n')

print("\ninput file created for run !!\n")

fid.close()
	

