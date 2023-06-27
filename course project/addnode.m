function addnode(tree,b1,b2,c1,k,gl,gr,psil,psir,tildef)
%将一个区间细化后的改变
m=size(tree.nodes,2);
n=size(tree.leaf,2);
i=1;
while ~isequal(tree.leaf(i).edge,[b1,b2])
    i=i+1;
end
l=tree.leaf(i).value;
new1=node(m+1,b1,c1);
new2=node(m+2,c1,b2);
new1.parent=tree.nodes(l);
new2.parent=tree.nodes(l);
new1.alpha(1)=CC(k,@(t) gl(t).*soldiscr(psil,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
new1.alpha(2)=CC(k,@(t) gr(t).*soldiscr(psil,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
new1.beta(1)=CC(k,@(t) gl(t).*soldiscr(psir,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
new1.beta(2)=CC(k,@(t) gr(t).*soldiscr(psir,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
new1.delta(1)=CC(k,@(t) gl(t).*soldiscr(tildef,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
new1.delta(2)=CC(k,@(t) gr(t).*soldiscr(tildef,k,t,b1,c1,gl,gr,psil,psir),b1,c1);
new2.alpha(1)=CC(k,@(t) gl(t).*soldiscr(psil,k,t,c1,b2,gl,gr,psil,psir),c1,b2);
new2.alpha(2)=CC(k,@(t) gr(t).*soldiscr(psil,k,t,c1,b2,gl,gr,psil,psir),c1,b2);
new2.beta(1)=CC(k,@(t) gl(t).*soldiscr(psir,k,t,c1,b2,gl,gr,psil,psir),c1,b2);
new2.beta(2)=CC(k,@(t) gr(t).*soldiscr(psir,k,t,c1,b2,gl,gr,psil,psir),c1,b2);
new2.delta(1)=CC(k,@(t) gl(t).*soldiscr(tildef,k,t,c1,b2,gl,gr,psil,psir),c1,b2);
new2.delta(2)=CC(k,@(t) gr(t).*soldiscr(tildef,k,t,c1,b2,gl,gr,psil,psir),c1,b2);
tree.nodes(l).left=new1;
tree.nodes(l).right=new2;
tree.nodes(l).alpha=[];
tree.nodes(l).beta=[];
tree.nodes(l).delta=[];
tree.nodes=[tree.nodes,new1,new2];
tree.leaf=[tree.leaf(1:i-1),new1,new2,tree.leaf(i+1:n)];