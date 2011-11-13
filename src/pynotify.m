function [status, result] = pynotify(title, msg)
    cmd = sprintf('python notifier.py "%s" "%s"', title, msg);
    [status, result] = system(cmd);
end
