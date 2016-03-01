classdef infos < handle
    properties
        data
    end
    
    methods
        function in = infos(inf)
            in.data=inf;
        end
    
        function varargout = subsref(in,S)
            varargout{:} = builtin('subsref', in.data, S);
%             varargout{:} = builtin('subsref', in, S);
        end

        function in = subsasgn(in,S,B)
%             in.data = builtin('subsasgn', in.data, S, B);
            in.data = builtin('subsasgn', in.data, S, B);
        end

        function data=getData(in)
            data=in.data;
        end
    end
end