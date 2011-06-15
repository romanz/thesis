function logger(level, msg, varargin)
if ~isempty(varargin)
    msg = sprintf(msg, varargin{:});
end
t = sprintf('%10.3f', toc);
fprintf('[%s] %s : %s\n', t, level, msg)
drawnow; % Update MATLAB GUI
