
L = 20

val_arr = rep(0,L)

i_beg = 10000
i_end = 11000

for (i in i_beg:i_end){

  N = i
  D = c(5,1)
  
  eigenValues_with_N = D^2/(4*N*N)
  
  dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  dRet = 0.0
  hyper0F1_val = .C("hyper_2by2_R",cc=1.5,eigenValues_with_N,dRet)[[3]]

  j = i-i_beg+1
  val_arr[j] = N*log(hyper0F1_val)
}

plot(exp(val_arr),type='l')
