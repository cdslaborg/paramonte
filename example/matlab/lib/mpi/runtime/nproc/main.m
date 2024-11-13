cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.lib.mpi.runtime.nproc()')
disp( pm.lib.mpi.runtime.nproc() )

disp("")
disp('pm.lib.mpi.runtime.nproc([])')
disp( pm.lib.mpi.runtime.nproc([]) )

disp("")
disp('pm.lib.mpi.runtime.nproc("impi")')
disp( pm.lib.mpi.runtime.nproc("impi") )

disp("")
disp('pm.lib.mpi.runtime.nproc("INTEL")')
disp( pm.lib.mpi.runtime.nproc("INTEL") )

disp("")
disp('pm.lib.mpi.runtime.nproc("intel")')
disp( pm.lib.mpi.runtime.nproc("intel") )

disp("")
disp('pm.lib.mpi.runtime.nproc("ompi")')
disp( pm.lib.mpi.runtime.nproc("ompi") )

disp("")
disp('pm.lib.mpi.runtime.nproc("openmpi")')
disp( pm.lib.mpi.runtime.nproc("openmpi") )

disp("")
disp('pm.lib.mpi.runtime.nproc("open-mpi")')
disp( pm.lib.mpi.runtime.nproc("open-mpi") )

disp("")
disp('pm.lib.mpi.runtime.nproc("mpich")')
disp( pm.lib.mpi.runtime.nproc("mpich") )

disp("")
disp('pm.lib.mpi.runtime.nproc("unknown")')
disp( pm.lib.mpi.runtime.nproc("unknown") )