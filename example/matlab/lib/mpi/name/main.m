cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.lib.mpi.name("impi")')
disp( pm.lib.mpi.name("impi") )
disp("")

disp('pm.lib.mpi.name("INTEL")')
disp( pm.lib.mpi.name("INTEL") )
disp("")

disp('pm.lib.mpi.name("intel")')
disp( pm.lib.mpi.name("intel") )
disp("")

disp('pm.lib.mpi.name("ompi")')
disp( pm.lib.mpi.name("ompi") )
disp("")

disp('pm.lib.mpi.name("openmpi")')
disp( pm.lib.mpi.name("openmpi") )
disp("")

disp('pm.lib.mpi.name("open-mpi")')
disp( pm.lib.mpi.name("open-mpi") )
disp("")

disp('pm.lib.mpi.name("mpich")')
disp( pm.lib.mpi.name("mpich") )
disp("")

disp('pm.lib.mpi.name("unknown")')
disp( pm.lib.mpi.name("unknown") )