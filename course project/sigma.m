function w=sigma(tree,k,x,gl,gr,psil,psir,tildef)
%计算给定二叉树后，积分方程在每个叶子结点上的解
s=size(tree.leaf,2);
b=zeros(1,s+1);
for i=1:s
    b(i)=tree.leaf(i).edge(1);
end
b(s+1)=tree.leaf(s).edge(2);
j=1;
while b(j)<x
    j=j+1;
end
j=max(j,2);
w=tree.leaf(j-1).mu(1)*soldiscr(tildef,k,x,b(j-1),b(j),gl,gr,psil,psir)+tree.leaf(j-1).mu(2)*soldiscr(psil,k,x,b(j-1),b(j),gl,gr,psil,psir)+tree.leaf(j-1).mu(3 ...
    )*soldiscr(psir,k,x,b(j-1),b(j),gl,gr,psil,psir);