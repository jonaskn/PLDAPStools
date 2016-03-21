function varargout = cellfun(varargin)
warning('OFF','MATLAB:dispatcher:nameConflict')
    varargout = cell(1,max(1,nargout));
    if nargin<2 || iscell(varargin{2})
        [varargout{:}] = builtin('cellfun',varargin{:});
    elseif isa(varargin{2}, 'lazyload.lazyload')
%         get the data
        fun=varargin{1};
        dat=cell(nargout,length(varargin{2}));
        for iField=1:length(varargin{2})
            [dat{:,iField}] = fun(varargin{2}{iField});
        end
    end
warning('ON','MATLAB:dispatcher:nameConflict')
end