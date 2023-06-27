classdef binary<handle
    properties
        nodes
        leaf
    end
    methods
        function obj=binary(n1)
            obj.nodes=[n1];
            obj.leaf=[n1];
        end
    end
end