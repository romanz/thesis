function list = glob(pattern)
    [path] = fileparts(pattern);
    d = dir(pattern);
    list = cell(1, numel(d));
    for k = 1:numel(d)
        list{k} = fullfile(path, d(k).name);
    end
end