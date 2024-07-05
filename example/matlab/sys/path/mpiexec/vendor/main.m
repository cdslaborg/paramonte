cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.vendor()')
                pm.sys.path.mpiexec.vendor()

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.vendor(pm.sys.path.mpiexec.which())')
                pm.sys.path.mpiexec.vendor(pm.sys.path.mpiexec.which())

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.vendor("paramonte")')
                pm.sys.path.mpiexec.vendor("paramonte")
assert(pm.sys.path.mpiexec.vendor("paramonte") == "")