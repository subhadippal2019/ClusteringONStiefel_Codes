function val=density1(y,D,S,dimIndex,N,bLog)
    %mex mhg.c
    %mex logmhg.c
    %D, S
    %
    %dim2=setdiff([1,2],dimIndex)
    %(d1*S(1,1))-N*(mhg(30,2,[],1.5,d))^N
    
    val=zeros(length(y),1);
    for i = 1:length(y)
        x=y(i);
        alpha=1;
        beta=0;
        D(dimIndex)=x;
        eigenValues=D.^2/4;
        if(bLog==1)
            val(i)=(alpha-1)*log(x)+(x*(S(dimIndex)-beta))- N * logmhg(60,2,[],1.5,eigenValues);
        else
            val(i)=x^(alpha-1)*exp(x*(S(dimIndex)-beta))/(mhg(10,2,[],1.5,eigenValues))^N;
        end
    end 
end