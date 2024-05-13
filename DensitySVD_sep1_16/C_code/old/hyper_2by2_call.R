


  dyn.load("hyper_2by2_R.so")
  val=seq(0,length(y));
  for i in 1:length(y)
    x=y(i);
    alpha=1;
    beta=0;
    D(dimIndex)=x;
    eigenValues=D.^2/4;

    if(bLog==1){
      val(i)=(alpha-1)*log(x)+(x*(S(dimIndex)-beta))- N * logmhg(60,2,[],1.5,eigenValues);
    }else{}
      val(i)=x^(alpha-1)*exp(x*(S(dimIndex)-beta))/(mhg(10,2,[],1.5,eigenValues))^N;

 
}