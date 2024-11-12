cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.lib.mpi.nproc()')
disp( pm.lib.mpi.nproc() )

disp("")
disp('pm.lib.mpi.nproc([])')
disp( pm.lib.mpi.nproc([]) )

disp("")
disp('pm.lib.mpi.nproc("impi")')
disp( pm.lib.mpi.nproc("impi") )

disp("")
disp('pm.lib.mpi.nproc("INTEL")')
disp( pm.lib.mpi.nproc("INTEL") )

disp("")
disp('pm.lib.mpi.nproc("intel")')
disp( pm.lib.mpi.nproc("intel") )

disp("")
disp('pm.lib.mpi.nproc("ompi")')
disp( pm.lib.mpi.nproc("ompi") )

disp("")
disp('pm.lib.mpi.nproc("openmpi")')
disp( pm.lib.mpi.nproc("openmpi") )

disp("")
disp('pm.lib.mpi.nproc("open-mpi")')
disp( pm.lib.mpi.nproc("open-mpi") )

disp("")
disp('pm.lib.mpi.nproc("mpich")')
disp( pm.lib.mpi.nproc("mpich") )

disp("")
disp('pm.lib.mpi.nproc("unknown")')
disp( pm.lib.mpi.nproc("unknown") )