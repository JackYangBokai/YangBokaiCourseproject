classdef node<handle
    %定义结点类
    properties
        value
        edge
        left
        right
        parent
        alpha
        beta
        delta
        mu
    end

    methods
        function obj=node(m,b1,b2)
        %创建结点
        obj.value=m;
        obj.edge=[b1,b2];
        obj.left=[];
        obj.right=[];
        obj.parent=[];
        obj.alpha=[];
        obj.beta=[];
        obj.delta=[];
        obj.mu=[];
        end
    end
end