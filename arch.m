s = datestr(now, 'yyyymmddHHMMSS');
s = sprintf('archive --prefix electrokinetics/ -o %s.zip HEAD', s);
git(s);
