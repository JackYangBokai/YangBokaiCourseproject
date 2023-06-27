function w=wholeche(f,k,b,x)
%给定一个grid和对应Chebyshev nodes上的值后，对应的插值
i=1;
while x>b(i)
    i=i+1;
end
i=max(i,2);
w=cheint(f((i-2)*k+1:(i-1)*k),k,b(i-1),b(i),x);