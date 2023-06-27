function w=soldiscr(f,k,x,b1,b2,g1,g2,psi1,psi2) 
%计算在一个区间上算子P离散化后的逆
t=zeros(1,k);
for i=1:k
    t(i)=(b2-b1)/2*cos((2*k-2*i+1)/(2*k)*pi)+(b1+b2)/2;
end
A=zeros(k,k);
for j=1:k
    h1=@(t) g1(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
    h2=@(t) g2(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
    [w1,c1]=intl(h1,k,b1,b2,t(1));
    [w2,c2]=intr(h2,k,b1,b2,t(1));
    e1=zeros(1,k);
    e2=zeros(1,k);
    e1(1)=w1;
    e2(1)=w2;
    for i=2:k
        for l=1:k
            e1(i)=e1(i)+(b2-b1)/2*c1(l)*cos((l-1)*acos(2*(t(i)-b1)/(b2-b1)-1));
            e2(i)=e2(i)+(b2-b1)/2*c2(l)*cos((l-1)*acos(2*(t(i)-b1)/(b2-b1)-1));
        end
    end
    for i=1:k
        A(i,j)=cos((j-1)*acos(2*(t(i)-b1)/(b2-b1)-1))+psi1(t(i))*e1(i)+psi2(t(i))*e2(i);
    end
end
c=f(t);
a=linsolve(A,c');
w=0;
for i=1:k
    w=w+a(i)*cos((i-1)*acos(2*(x-b1)/(b2-b1)-1));
end
 

