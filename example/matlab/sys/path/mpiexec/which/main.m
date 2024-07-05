cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.which()')
                pm.sys.path.mpiexec.which()

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.which("intel")')
                pm.sys.path.mpiexec.which("intel")

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.which("mpich")')
                pm.sys.path.mpiexec.which("mpich")

pm.matlab.show()
pm.matlab.show('pm.sys.path.mpiexec.which("openmpi")')
                pm.sys.path.mpiexec.which("openmpi")