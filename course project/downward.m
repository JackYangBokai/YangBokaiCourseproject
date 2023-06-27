function downward(tree,node)
%给定一颗二叉树，计算coupling coefficient的向下传递
if isempty(node.left)
    return;
end
l=node.value;
[tree.nodes(l).left.mu,tree.nodes(l).right.mu]=coupcoe(tree.nodes(l).mu, ...
    tree.nodes(l).left.alpha,tree.nodes(l).left.beta,tree.nodes(l).left.delta, ...
    tree.nodes(l).right.alpha,tree.nodes(l).right.beta,tree.nodes(l).right.delta);
downward(tree,node.left);
downward(tree,node.right);