cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.lib.mpi.verify()')
disp( pm.lib.mpi.verify() )

disp("")
disp('pm.lib.mpi.verify("Intel")')
disp( pm.lib.mpi.verify("Intel") )

disp("")
disp('pm.lib.mpi.verify("MPICH")')
disp( pm.lib.mpi.verify("MPICH") )

disp("")
disp('pm.lib.mpi.verify("OpenMPI")')
disp( pm.lib.mpi.verify("OpenMPI") )

disp("")
disp('pm.lib.mpi.verify("all")')
disp( pm.lib.mpi.verify("all") )

disp("")
disp('pm.lib.mpi.verify("any")')
disp( pm.lib.mpi.verify("any") )