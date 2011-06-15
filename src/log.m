function log(level, msg, varargin)
if ~isempty(varargin)
    msg = sprintf(msg, varargin{:});
end
fprintf('%-10s : %s\n', level, msg)