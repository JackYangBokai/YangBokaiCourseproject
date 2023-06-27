function [d,e]=coupcoe(b,alphad,betad,deltad,alphae,betae,deltae)
%计算coupling coefficient的向下传递
d=zeros(1,3);
e=zeros(1,3);
d(1)=b(1);
e(1)=b(1);
d(2)=b(2);
e(3)=b(3);
A=[1,alphae(2);betad(1),1];
f=[b(3)*(1-betae(2))-b(1)*deltae(2);b(2)*(1-alphad(1))-b(1)*deltad(1)];
g=linsolve(A,f);
d(3)=g(1);
e(2)=g(2);