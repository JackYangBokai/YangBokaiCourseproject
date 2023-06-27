function [w1,w2]=program3(f,f1,f2,a,b,c,k,zetal1,zetal2,zetar1,zetar2,u0,u01,u02,eps,con,TOL)
%实现第三个程序，f1为f关于u'的导数，f2为f关于u的导数
%边值条件使用前两问中齐次的形式,u0为解的initial guess，u01,u02为解的一阶导与二阶导的initial guess,con代替program2中的C

%在这里，我们没有运用第二篇论文的做法，而是仿照将第二篇论文的想法构造了类似的residual mapping，并且得到
%迭代方程为v''-f1*v-f2=-(u_n''-f(u_n',u_n,x)),这里u_n是第n步的解，而v=u_{n+1}-u_n为差
test=1;
u=u0;
u1=u01;
u2=u02;
while test>eps
    v=u;
    v1=u1;
    v2=u2;
    p=@(t) -f1(v1(t),v(t),t);
    q=@(t) -f2(v1(t),v(t),t);
    g=@(t) -v2(t)+f(v1(t),v(t),t);
    [del,del1,tree]=program2(p,q,g,a,b,c,k,zetal1,zetal2,zetar1,zetar2,0,0,con,TOL);
    %将计算结果储存在Chebyshev结点上再进行插值，以防每次后面计算再次进行program2的计算过程
    si=size(tree.leaf,2);
    b1=zeros(1,si+1);
    for i=1:si
        b1(i)=tree.leaf(i).edge(1);
    end
    b1(si+1)=tree.leaf(si).edge(2);
    t=zeros(1,si*k);
    deld=zeros(1,si*k);
    del1d=zeros(1,si*k);
    for i=1:si
        for j=1:k
            t((i-1)*k+j)=(b1(i+1)-b1(i))/2*cos((2*k-2*j+1)/(2*k)*pi)+(b1(i+1)+b1(i))/2;
            deld((i-1)*k+j)=del(t((i-1)*k+j));
            del1d((i-1)*k+j)=del1(t((i-1)*k+j));
        end
    end
    del=@(t) wholeche(deld',k,b1,t);
    del1=@(t) wholeche(del1d',k,b1,t);
    %计算新的结果
    u=@(t) v(t)+del(t);
    u1=@(t) v1(t)+del1(t);
    u2=@(t) v2(t)+f1(v1(t),v(t),t).*del1(t)+f2(v1(t),v(t),t).*del(t)+g(t);
    %计算test
    e1=integral(@(t) del(t).^2,a,c);
    e2=integral(@(t) v(t).^2,a,c,"ArrayValued",true);
    test=sqrt(e1/e2);
end
w1=u;
w2=u1;