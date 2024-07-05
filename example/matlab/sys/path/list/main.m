cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.sys.path.list()')
pm.matlab.show( pm.sys.path.list() )

pm.matlab.show()
pm.matlab.show('pm.sys.path.list("*.out*")')
pm.matlab.show( pm.sys.path.list("*.out*") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.list("..")')
pm.matlab.show( pm.sys.path.list("..") )