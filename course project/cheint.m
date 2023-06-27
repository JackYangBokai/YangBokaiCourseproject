function [w,a]=cheint(f,k,b1,b2,x)
%计算f的切比雪夫插值结果,其中a是切比雪夫多项式前的系数
t=zeros(1,k);
for i=1:k
    t(i)=cos((2*k-2*i+1)/(2*k)*pi);
end
a=zeros(1,k);
for i=2:k
    a(i)=2/k*cos((i-1)*acos(t))*f;
end
a(1)=1/k*sum(f);
w=0;
for i=1:k
    w=w+a(i)*cos((i-1)*acos(2*(x-b1)/(b2-b1)-1));
end
