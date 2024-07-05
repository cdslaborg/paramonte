cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.sys.path.glob("*.m")')
pm.matlab.show( pm.sys.path.glob("*.m") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.glob("example/array/**")')
pm.matlab.show( pm.sys.path.glob("example/array/**") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.glob("example/array/**.m")')
pm.matlab.show( pm.sys.path.glob("example/array/**.m") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.glob("example/array/*.m")')
pm.matlab.show( pm.sys.path.glob("example/array/*.m") )

pm.matlab.show()
pm.matlab.show('pm.sys.path.glob(".*")')
pm.matlab.show( pm.sys.path.glob(".*") )

pm.matlab.show()
pm.matlab.show('[paths, isdir] = pm.sys.path.glob("**"); paths(~isdir) % get all files in directory tree.')
                [paths, isdir] = pm.sys.path.glob("**"); paths(~isdir) % get all files in directory tree.