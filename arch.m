s = datestr(now, 'yyyymmddHHMM');
s = sprintf('archive -o %s.zip HEAD', s);
git(s);
