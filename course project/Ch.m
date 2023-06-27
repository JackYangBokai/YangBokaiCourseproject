function w=Ch(k,b1,b2)
%写出chebyshev integration的矩阵的具体形式,其中该矩阵作用在函数在切比雪夫结点上的值
w=zeros(1,k);
t=zeros(1,k);
for i=1:k
    t(i)=cos((2*k-2*i+1)/(2*k)*pi);
end
for i=1:k
    a=zeros(k,1);
    a(i)=1;
    [~,f]=cheint(a,k,b1,b2,t);
    w(i)=0;
    for l=1:2:k
        w(i)=w(i)+(b2-b1)/2*f(l)*2/(1-(l-1)^2);
    end
end