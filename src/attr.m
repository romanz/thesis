function vals = attr(objs, name)
    vals = cell(size(objs));
    for k = 1:numel(vals)
        vals{k} = objs{k}.(name);
    end
    vals = cell2mat(vals);
end