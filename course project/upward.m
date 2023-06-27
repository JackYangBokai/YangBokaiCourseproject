function upward(tree,node)
%给定一颗二叉树，计算内积的向上传递
if isempty(node.left)
    return;
end
upward(tree,node.left);
upward(tree,node.right);
l=node.value;
[tree.nodes(l).alpha,tree.nodes(l).beta,tree.nodes(l).delta]=innpro( ...
    tree.nodes(l).left.alpha,tree.nodes(l).right.alpha,tree.nodes(l).left.beta, ...
    tree.nodes(l).right.beta,tree.nodes(l).left.delta,tree.nodes(l).right.delta);