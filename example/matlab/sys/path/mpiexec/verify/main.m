cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.sys.path.mpiexec.verify()')
disp( pm.sys.path.mpiexec.verify() )

disp("")
disp('pm.sys.path.mpiexec.verify("Intel")')
disp( pm.sys.path.mpiexec.verify("Intel") )

disp("")
disp('pm.sys.path.mpiexec.verify("MPICH")')
disp( pm.sys.path.mpiexec.verify("MPICH") )

disp("")
disp('pm.sys.path.mpiexec.verify("OpenMPI")')
disp( pm.sys.path.mpiexec.verify("OpenMPI") )

disp("")
disp('pm.sys.path.mpiexec.verify("all")')
disp( pm.sys.path.mpiexec.verify("all") )

disp("")
disp('pm.sys.path.mpiexec.verify("any")')
disp( pm.sys.path.mpiexec.verify("any") )