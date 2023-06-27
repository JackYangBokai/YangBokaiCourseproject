function [alphab,betab,deltab]=innpro(alphad,alphae,betad,betae,deltad,deltae)
%计算内积系数的向上传递
del=1-alphae(2)*betad(1);
alphab=zeros(1,2);
betab=zeros(1,2);
deltab=zeros(1,2);
alphab(1)=(1-alphae(1))*(alphad(1)-betad(1)*alphae(2))/del+alphae(1);
alphab(2)=alphae(2)*(1-betad(2))*(1-alphad(1))/del+alphad(2);
betab(1)=betad(1)*(1-betae(2))*(1-alphae(1))/del+betae(1);
betab(2)=(1-betad(2))*(betae(2)-betad(1)*alphae(2))/del+betad(2);
deltab(1)=(1-alphae(1))*(deltad(1)-betad(1)*deltae(2))/del+deltae(1);
deltab(2)=(1-betad(2))*(deltae(2)-alphae(2)*deltad(1))/del+deltad(2);