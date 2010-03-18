python -c "import os, time; os.system('hg arch -t zip numerics_%%s.zip' %% time.strftime('%%Y.%%m.%%d_%%H.%%M.%%S'))"
pause