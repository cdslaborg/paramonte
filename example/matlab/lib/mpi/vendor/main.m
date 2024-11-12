cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.lib.mpi.vendor("impi")')
disp( pm.lib.mpi.vendor("impi") )
disp("")

disp('pm.lib.mpi.vendor("INTEL")')
disp( pm.lib.mpi.vendor("INTEL") )
disp("")

disp('pm.lib.mpi.vendor("intel")')
disp( pm.lib.mpi.vendor("intel") )
disp("")

disp('pm.lib.mpi.vendor("ompi")')
disp( pm.lib.mpi.vendor("ompi") )
disp("")

disp('pm.lib.mpi.vendor("openmpi")')
disp( pm.lib.mpi.vendor("openmpi") )
disp("")

disp('pm.lib.mpi.vendor("open-mpi")')
disp( pm.lib.mpi.vendor("open-mpi") )
disp("")

disp('pm.lib.mpi.vendor("mpich")')
disp( pm.lib.mpi.vendor("mpich") )
disp("")

disp('pm.lib.mpi.vendor("unknown")')
disp( pm.lib.mpi.vendor("unknown") )