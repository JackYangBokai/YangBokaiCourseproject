function w=monitor(f,k,b1,b2,g1,g2,psi1,psi2)
%计算误差估计需要的monitor function
t=zeros(1,k);
for i=1:k
    t(i)=(b2-b1)/2*cos((2*k-2*i+1)/(2*k)*pi)+(b1+b2)/2;
end
A=zeros(k,k);
for i=1:k
    for j=1:k
        h1=@(t) g1(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
        h2=@(t) g2(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
        A(i,j)=cos((j-1)*acos(2*(t(i)-b1)/(b2-b1)-1))+psi1(t(i))*intl(h1, ...
            k,b1,b2,t(i))+psi2(t(i))*intr(h2,k,b1,b2,t(i));
    end
end
c=f(t);
a=linsolve(A,c');
w=abs(a(k-1))+abs(a(k)-a(k-2));