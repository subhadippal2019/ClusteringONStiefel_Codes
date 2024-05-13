bivariate_NR_method_d1_d2 <-function(p,g1,g2,max_iter){

 ### input D is a initial value from Mardia's book for small value estimate
 ### d1 = p*g_1; d2 = p*g2  
 
 print(paste0("inside bivariate NR method for updating d1 and d2"))
 
 #dyn.load("./C_code/functions_hyper_2by2_R.so") 

 nrow_D = p 
 ncol_D = p
 ### compute H
 H = matrix(rep(0.0,nrow_D*ncol_D),ncol=ncol_D)
 
 G = rep(0.0,2)
 G[1] = g1
 G[2] = g2
 
 ### D = p*G
 D = rep(0.0,2)
 D_old = rep(0.0,2)
 D_new = rep(0.0,2)
 
 D[1] = p*g1
 D[2] = p*g2
 
 D_old = D
 D_new = D_old
 
 for (i in  1:max_iter){
   
   D_old = D_new
  
   eigenValues=D_old^2/4
   dRet = 0.0
   
   hyper0F1_val = .C("hyper_2by2_R",cc=1.5,eigenValues,dRet)[[3]]
   
   Y1 = .C("partial_d1_hyper_2by2_R",1.5,eigenValues,dRet)[[3]]
   YY1 = .C("partial_d1_partial_d1_hyper_2by2_R",1.5,eigenValues,dRet)[[3]]
   
   YY12 = .C("partial_d1_partial_d2_hyper_2by2_R",1.5,eigenValues,dRet)[[3]]
   YY21 = YY12
   
   Y2 = .C("partial_d2_hyper_2by2_R",1.5,eigenValues,dRet)[[3]]
   YY2 = .C("partial_d2_partial_d2_hyper_2by2_R",1.5,eigenValues,dRet)[[3]]
   
   H[1,1] = YY1/hyper0F1_val - (Y1/hyper0F1_val)*(Y1/hyper0F1_val)
   H[1,2] = YY12/hyper0F1_val - (Y1/hyper0F1_val)*(Y2/hyper0F1_val)
   H[2,1] = H[1,2]
   H[2,2] = YY2/hyper0F1_val - (Y2/hyper0F1_val)*(Y2/hyper0F1_val)
   
   F = rep(0.0,2)
   F[1] = Y1/hyper0F1_val
   F[2] = Y2/hyper0F1_val
   
   H_inv = solve(H)
   
   D_new = D_old - H_inv%*%(F-G)
   
   ##print(D_new)
 }
 
 return(as.vector(D_new)) 
}