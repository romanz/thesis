function conf = defaults(conf, varargin)
    if nargin < 1 || isempty(conf)
        conf = {};
    end
    if iscell(conf)
        conf = struct(conf{:});
    end

    args = struct(varargin{:});
    names = fieldnames(args);
    for k = 1:numel(names)
        if ~isfield(conf, names{k}) % Copy only missing fields
            conf.(names{k}) = args.(names{k});
        end
    end
    
end
