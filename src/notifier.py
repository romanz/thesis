import pynotify
import sys

title, msg = sys.argv[1:]
assert( pynotify.init('MATLAB') ) 
n = pynotify.Notification(title, msg)
n.show()
