function [w1,w2,tree]=program2(p,q,f,a,b,c,k,zetal1,zetal2,zetar1,zetar2,gammal,gammar,C,TOL)
%实现第二个程序（full algorithm），其中输出w1为原方程的解，w2为该解的导数
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
%构造一开始的二叉树结构
m=size(b,2);
n1=node(1,b(1),b(m));
n1.mu=[1,0,0];
tree=binary(n1);
for i=1:m-2
    addnode(tree,b(1),b(m-i+1),b(m-i),k,gl,gr,psil,psir,tildef);
end
%进行系数的计算
upward(tree,tree.nodes(1));
downward(tree,tree.nodes(1));
%进行mesh refinement
test=1;
while test>TOL
    oldtree.nodes=tree.nodes;
    oldtree.leaf=tree.leaf;
    v=@(t) origin(oldtree,ui,k,gl,gr,psil,psir,tildef,s,t);
    si=size(tree.leaf,2);
    mo=zeros(1,si);
    val=zeros(1,si);
    con=zeros(1,si);
    for i=1:si
        app=@(t) tree.leaf(i).mu(1)*tildef(t)+tree.leaf(i).mu(2)*psil(t ...
            )+tree.leaf(i).mu(3)*psir(t);
        mo(i)=monitor(app,k,tree.leaf(i).edge(1),tree.leaf(i).edge(2),gl,gr,psil,psir);
        val(i)=tree.leaf(i).value;
    end
    sdiv=max(mo)/(2^C);
    for i=1:si
        if i>1 && con(i-1)>0
            continue;
            %如果对前一个结点使用了merge操作，直接跳过该结点
        elseif mo(i)>=sdiv
            addnode(tree,tree.nodes(val(i)).edge(1),tree.nodes(val(i)).edge(2), ...
                1/2*(tree.nodes(val(i)).edge(1)+tree.nodes(val(i)).edge(2)), ...
                k,gl,gr,psil,psir,tildef);
        elseif i==si
            continue;
            %防止对最后一个结点进行merge操作
        elseif mo(i)+mo(i+1)<sdiv/(2^k) && tree.nodes(val(i)).parent==tree.nodes(val(i+1)).parent
            merge(tree,tree.nodes(val(i)),tree.nodes(val(i+1)), ...
                k,gl,gr,psil,psir,tildef);
            con(i)=1;
            for j=i+2:si
                if val(j)>val(i)+1
                    val(j)=val(j)-2;
                end
            end
        end
    end
    upward(tree,tree.nodes(1));
    downward(tree,tree.nodes(1));
    %计算test，以决定继续进行mesh refinement还是对所有区间进行最后的二分
    si=size(tree.leaf,2);
    e1=0;
    e2=0;
    u=@(t) origin(tree,ui,k,gl,gr,psil,psir,tildef,s,t);
    dif=@(t) (u(t)-v(t)).^2;
    sum=@(t) (u(t)+v(t)).^2;
    for i=1:si
        e1=e1+integral(dif,tree.leaf(i).edge(1),tree.leaf(i).edge(2));
        e2=e2+integral(sum,tree.leaf(i).edge(1),tree.leaf(i).edge(2));
    end
    test=sqrt(e1/e2);
end
%最后一步：对每个叶子结点再进行划分
si=size(tree.leaf,2);
val=zeros(1,si);
for i=1:si
    val(i)=tree.leaf(i).value;
end
for i=1:si
    addnode(tree,tree.nodes(val(i)).edge(1),tree.nodes(val(i)).edge(2), ...
        1/2*(tree.nodes(val(i)).edge(1)+tree.nodes(val(i)).edge(2)), ...
        k,gl,gr,psil,psir,tildef);
end
upward(tree,tree.nodes(1));
downward(tree,tree.nodes(1));
%计算结果
w1=@(t) origin(tree,ui,k,gl,gr,psil,psir,tildef,s,t);
w2=@(t) origin1(tree,ui1,k,gl1,gr1,gl,gr,psil,psir,tildef,s,t);