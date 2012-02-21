function s = tolatex(s)
s = simple(s);
s = latex(s);
s = strrep(s, '\mathrm{U', '{\cU_');
s = strrep(s, '\, a', '\alpha');
s = strrep(s, '\,', '');
s = strrep(s, '\ln', '\log');
s = strrep(s, '{g}', '{\gamma}');
s = strrep(s, '\!', '');
s = strrep(s, '\left(t\right)}^3', '^3\theta}');
s = strrep(s, '\left(t\right)}^2', '^2\theta}');
s = strrep(s, '\left(t\right)', '\theta');