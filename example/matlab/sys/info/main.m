cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.sys.info()')
pm.matlab.show( pm.sys.info() )

pm.matlab.show()
pm.matlab.show('[str, cache] = pm.sys.info();')
                [str, cache] = pm.sys.info();
pm.matlab.show('cache')
pm.matlab.show( cache )
pm.matlab.show('str')
pm.matlab.show( str )