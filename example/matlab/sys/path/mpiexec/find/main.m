cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.find()')
                pm.sys.path.mpiexec.find()

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.find("intel")')
                pm.sys.path.mpiexec.find("intel")

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.find("mpich")')
                pm.sys.path.mpiexec.find("mpich")

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.find("openmpi")')
                pm.sys.path.mpiexec.find("openmpi")
