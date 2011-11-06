function v = dukhin_derjaguin(phi, g)
    v = 4 * log((exp((-phi-log(g))/2) + 1)/2) * diff(phi, 't');
end
