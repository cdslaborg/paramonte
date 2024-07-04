cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.lib.version()')
pm.matlab.show( pm.lib.version() )

pm.matlab.show()
pm.matlab.show('pm.lib.version("matlab")')
pm.matlab.show( pm.lib.version("matlab") )

pm.matlab.show()
pm.matlab.show('pm.lib.version("generic")')
pm.matlab.show( pm.lib.version("generic") )

pm.matlab.show()

pm.matlab.show()
pm.matlab.show('pm.lib.version([], "major")')
pm.matlab.show( pm.lib.version([], "major") )

pm.matlab.show()
pm.matlab.show('pm.lib.version([], "minor")')
pm.matlab.show( pm.lib.version([], "minor") )

pm.matlab.show()
pm.matlab.show('pm.lib.version([], "patch")')
pm.matlab.show( pm.lib.version([], "patch") )

pm.matlab.show()

pm.matlab.show()
pm.matlab.show('pm.lib.version("matlab", "major")')
pm.matlab.show( pm.lib.version("matlab", "major") )

pm.matlab.show()
pm.matlab.show('pm.lib.version("matlab", "minor")')
pm.matlab.show( pm.lib.version("matlab", "minor") )

pm.matlab.show()
pm.matlab.show('pm.lib.version("matlab", "patch")')
pm.matlab.show( pm.lib.version("matlab", "patch") )

pm.matlab.show()

pm.matlab.show()
pm.matlab.show('pm.lib.version("generic", "major")')
pm.matlab.show( pm.lib.version("generic", "major") )

pm.matlab.show()
pm.matlab.show('pm.lib.version("generic", "minor")')
pm.matlab.show( pm.lib.version("generic", "minor") )

pm.matlab.show()
pm.matlab.show('pm.lib.version("generic", "patch")')
pm.matlab.show( pm.lib.version("generic", "patch") )