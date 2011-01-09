function handle = progress(handle, p, msg)

if isempty(handle)
    handle.start = now;
    handle.fig = waitbar(0, '');
    handle.msg = msg;
    msg = '';
    handle.p = 0;
elseif isempty(p)
    close(handle.fig);
    return
% elseif p < min(1, handle.p + 1e-2)
%     return
end

if nargin >= 3 && ~isempty(msg)
    msg = sprintf('\n%s', msg);
end

handle.p = p;

runtime = datestr(now - handle.start, 'HH:MM:SS');
runtime = ['Runtime: ' runtime ' [seconds]'];

waitbar(handle.p, handle.fig, [runtime msg]);
s = [sprintf('%.0f%%', 100*handle.p) ' - ' handle.msg];
set(handle.fig, 'Name', s);
