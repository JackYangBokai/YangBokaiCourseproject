function w=origin(tree,ui,k,gl,gr,psil,psir,tildef,s,x)
%给定一棵二叉树，计算原方程的解
si=size(tree.leaf,2);
b=zeros(1,si+1);
for i=1:si
    b(i)=tree.leaf(i).edge(1);
end
b(si+1)=tree.leaf(si).edge(2);
j=1;
while b(j)<x
    j=j+1;
end
j=max(j,2);
[Jl,Jr]=precom(tree);
sig=@(t) sigma(tree,k,t,gl,gr,psil,psir,tildef);
w=ui(x)+gr(x)/s.*(Jl(j-1)+intl(@(t) gl(t).*sig(t),k,b(j-1),b(j),x))+gl(x)/s.*( ...
    Jr(j)+intr(@(t) gr(t).*sig(t),k,b(j-1),b(j),x));