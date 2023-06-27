function [Jl,Jr]=precom(tree)
%给定一颗二叉树后，类似论文3.A步骤提前计算定积分
s=size(tree.leaf,2);
Jl=zeros(1,s+1);
Jr=zeros(1,s+1);
Jl(1)=0;
Jr(s+1)=0;
for i=1:s
    Jl(i+1)=Jl(i)+tree.leaf(i).delta(1)*tree.leaf(i).mu(1)+tree.leaf(i).alpha(1 ...
        )*tree.leaf(i).mu(2)+tree.leaf(i).beta(1)*tree.leaf(i).mu(3);
end
for i=s+1:-1:2
    Jr(i-1)=Jr(i)+tree.leaf(i-1).delta(2)*tree.leaf(i-1).mu(1)+tree.leaf(i-1).alpha(2 ...
        )*tree.leaf(i-1).mu(2)+tree.leaf(i-1).beta(2)*tree.leaf(i-1).mu(3);
end