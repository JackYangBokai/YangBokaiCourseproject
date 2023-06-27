function u=program1(p,q,f,a,b,c,k,zetal1,zetal2,zetar1,zetar2,gammal,gammar)
%实现第一个程序（direct solver）,其中输出u第一行为原方程在Chebyshev结点上的解，第二行为该解的导数在Chebyshev结点上的解
%化为积分方程
lam=[zetal1*a+zetal2, zetal1;zetar1*c+zetar2, zetar1];
gam=[gammal;gammar];
l=linsolve(lam,gam);
ui=@(t) l(1)*t+l(2);
ui1=@(t) l(1);
tildef=@(t) f(t)-(p(t)*l(1)+q(t).*ui(t));
if abs(zetal1)>=abs(zetal2) || abs(zetar1)>=abs(zetar2)
    gl=@(t) zetal1*(t-a)-zetal2;
    gr=@(t) zetar1*(t-c)-zetar2;
    gl1=@(t) zetal1;
    gr1=@(t) zetar1;
    s=(zetar1*c+zetar2)*zetal1-(zetal1*a+zetal2)*zetar1;
else
    gl=@(t) zetal2*cosh(t-a)-zetal1*sinh(t-a);
    gr=@(t) zetar2*cosh(t-c)-zetar1*sinh(t-c);
    gl1=@(t) zetal2*sinh(t-a)-zetal1*cosh(t-a);
    gr1=@(t) zetar2*sinh(t-c)-zetar1*cosh(t-c);
    s=gl(a)*gr1(a)-gr(a)*gl1(a);
    q=@(t) q(t)+1;
end
psil=@(t) (p(t).*gr1(t)+q(t).*gr(t))/s;
psir=@(t) (p(t).*gl1(t)+q(t).*gl(t))/s;
%列出Chebyshev结点
n=length(b)-1;
t=zeros(1,n*k);
for i=1:n
    for j=1:k
         t((i-1)*k+j)=(b(i+1)-b(i))/2*cos((2*k-2*j+1)/(2*k)*pi)+(b(i+1)+b(i))/2;
    end
end
%计算具体矩阵
P=zeros(n*k,n*k);
for i=1:n
    for j=1:k
        T=(i-1)*k+j;
        for v=1:i-1
            for w=1:k
                e=zeros(k,1);
                e(w)=1;
                z=(v-1)*k+w;
                P(T,z)=psil(t(T))*gl(t(z))*Ch(k,b(v),b(v+1))*e;
            end
        end
        for v=i+1:n
            for w=1:k
                e=zeros(k,1);
                e(w)=1;
                z=(v-1)*k+w;
                P(T,z)=psir(t(T))*gr(t(z))*Ch(k,b(v),b(v+1))*e;
            end
        end
        for w=1:k
            z=(i-1)*k+w;
            A1=Isl(k,b(i),b(i+1));
            A2=Isr(k,b(i),b(i+1));
            P(T,z)=psil(t(T))*gl(t(z))*A1(j,w)+psir(t(T))*gr(t(z))*A2(j,w);
        end
    end
end
P=eye(n*k)+P;
a=tildef(t);
si=linsolve(P,a');
%已知sigma在chebyshev结点上的值，求原解
%其中利用第二个程序预先计算定积分的思想
Jl=zeros(1,n+1);
Jr=zeros(1,n+1);
Jl(1)=0;
Jr(n+1)=0;
for i=1:n
    Jl(i+1)=Jl(i)+Ch(k,b(i),b(i+1))*(gl(t((i-1)*k+1:i*k))'.*si((i-1)*k+1:i*k));
end
for i=n+1:-1:2
    Jr(i-1)=Jr(i)+Ch(k,b(i-1),b(i))*(gr(t((i-2)*k+1:(i-1)*k))'.*si((i-2)*k+1:(i-1)*k));
end
u=zeros(2,n*k);
for i=1:n
    sig=si((i-1)*k+1:i*k);
    g1=gl(t((i-1)*k+1:i*k))';
    g2=gr(t((i-1)*k+1:i*k))';
    for j=1:k
        e=zeros(1,k);
        e(j)=1;
        T=(i-1)*k+j;
        u(1,T)=ui(t(T))+gr(t(T))/s*(Jl(i)+e*Isl(k,b(i),b(i+1))*(g1.*sig))+gl(t(T))/s*(Jr(i+1)+e*Isr(k,b(i),b(i+1))*(g2.*sig));
        u(2,T)=ui1(t(T))+gr1(t(T))/s*(Jl(i)+e*Isl(k,b(i),b(i+1))*(g1.*sig))+gl1(t(T))/s*(Jr(i+1)+e*Isr(k,b(i),b(i+1))*(g2.*sig));
    end
end
