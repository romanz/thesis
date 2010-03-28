python -c "import os, time; os.system('hg arch -X *.cmd -t zip numerics_%%s.zip' %% time.strftime('%%Y.%%m.%%d'))"
pause