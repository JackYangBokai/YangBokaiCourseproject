function w=CC(k,f,b1,b2)
%用Chebyshev结点计算积分
t=zeros(1,k);
for i=1:k
    t(i)=cos((2*k-2*i+1)/(2*k)*pi);
end
a=zeros(1,k);
for i=2:k
    a(i)=2/k*sum(f((b2-b1)/2*t+(b2+b1)/2).*cos((i-1)*acos(t)));
end
a(1)=1/k*sum(f((b2-b1)/2*t+(b2+b1)/2));
w=0;
for i=1:2:k
    w=w+(b2-b1)/2*a(i)*2/(1-(i-1)^2);
end

