classdef State < handle
    methods(Abstract)
        handleRequest(obj, context);
    end
end