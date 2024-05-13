chk_timimg <-function(M){
  library(plyr)
  Rprof(interval=0.02)
  dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  t1=Sys.time()
  for(i in 1:M){
    dRet = 0.0
    D = c(21+rnorm(1),17+rnorm(1))
    eigenValues = D^2/4
    hyper0F1_val = .C("hyper_2by2_R",a=1.5,eigenValues,dRet)[[3]]
  }
  Sys.time()-t1
  Rprof(NULL)
  out=summaryRprof()
  return(out)
}

