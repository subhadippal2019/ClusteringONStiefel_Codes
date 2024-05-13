function S= hyperGeometric0F1_3by2and2times2_1(D,K,N)
if(isvector(D))
    x=D(1);
    y=D(2);
end 
S=0;
for k = 0:K
    for n =0:N
        S=S+((x*y)^k*(x+y)^n)/(factorial(k)*factorial(n)*pochhammer(.5,k)*pochhammer(1.5,2*k)*pochhammer(1.5+2*k,n));
    end
end
