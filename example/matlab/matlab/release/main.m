cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.matlab.release()')
pm.matlab.show( pm.matlab.release() )

pm.matlab.show()
pm.matlab.show('str2double(pm.matlab.release("year"))')
pm.matlab.show( str2double(pm.matlab.release("year")) )

pm.matlab.show()
pm.matlab.show('pm.matlab.release("minor")')
pm.matlab.show( pm.matlab.release("minor") )

pm.matlab.show()
pm.matlab.show('pm.matlab.release("season")')
pm.matlab.show( pm.matlab.release("season") )
