function merge(tree,node1,node2,k,gl,gr,psil,psir,tildef)
%消除两个parent相同的叶子节点
if ~isempty(node1.left)
    return;
end
if ~isempty(node2.left)
    return;
end
if node1.parent~=node2.parent
    return;
end
l1=node1.value;
l2=node2.value;
l=node1.parent.value;
i=1;
while ~isequal(tree.leaf(i).value,l1)
    i=i+1;
end
if l1>l2
    i=i-1;
end
n=size(tree.leaf,2);
b1=tree.nodes(l).edge(1);
c1=tree.nodes(l).edge(2);
tree.nodes(l).left=[];
tree.nodes(l).right=[];
tree.nodes(l).alpha(1)=CC(k,@(t) gl(t).*soldiscr(psil,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
tree.nodes(l).alpha(2)=CC(k,@(t) gr(t).*soldiscr(psil,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
tree.nodes(l).beta(1)=CC(k,@(t) gl(t).*soldiscr(psir,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
tree.nodes(l).beta(2)=CC(k,@(t) gr(t).*soldiscr(psir,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
tree.nodes(l).delta(1)=CC(k,@(t) gl(t).*soldiscr(tildef,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
tree.nodes(l).delta(2)=CC(k,@(t) gr(t).*soldiscr(tildef,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
tree.leaf=[tree.leaf(1:i-1),tree.nodes(l),tree.leaf(i+2:n)];
tree.nodes(l1)=[];
if l1>l2
    tree.nodes(l2)=[];
else
    tree.nodes(l2-1)=[];
end
si=size(tree.nodes,2);
for k=1:si
    tree.nodes(k).value=k;
end

