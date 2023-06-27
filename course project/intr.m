function [w,c]=intr(f,k,b1,b2,x) 
%实现一般区间上的right integration operator
t=zeros(1,k);
for i=1:k
    t(i)=cos((2*k-2*i+1)/(2*k)*pi);
end
a=zeros(1,k);
for i=2:k
    a(i)=2/k*sum(f((b2-b1)/2*t+(b2+b1)/2).*cos((i-1)*acos(t)));
end
a(1)=1/k*sum(f((b2-b1)/2*t+(b2+b1)/2));
c=intsr(a);
w=0;
for i=1:k
    w=w+(b2-b1)/2*c(i)*cos((i-1)*acos(2*(x-b1)/(b2-b1)-1));
end