cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.lib.mpi.runtime.rankp1()')
disp( pm.lib.mpi.runtime.rankp1() )

disp("")
disp('pm.lib.mpi.runtime.rankp1([])')
disp( pm.lib.mpi.runtime.rankp1([]) )

disp("")
disp('pm.lib.mpi.runtime.rankp1("impi")')
disp( pm.lib.mpi.runtime.rankp1("impi") )

disp("")
disp('pm.lib.mpi.runtime.rankp1("INTEL")')
disp( pm.lib.mpi.runtime.rankp1("INTEL") )

disp("")
disp('pm.lib.mpi.runtime.rankp1("intel")')
disp( pm.lib.mpi.runtime.rankp1("intel") )

disp("")
disp('pm.lib.mpi.runtime.rankp1("ompi")')
disp( pm.lib.mpi.runtime.rankp1("ompi") )

disp("")
disp('pm.lib.mpi.runtime.rankp1("openmpi")')
disp( pm.lib.mpi.runtime.rankp1("openmpi") )

disp("")
disp('pm.lib.mpi.runtime.rankp1("open-mpi")')
disp( pm.lib.mpi.runtime.rankp1("open-mpi") )

disp("")
disp('pm.lib.mpi.runtime.rankp1("mpich")')
disp( pm.lib.mpi.runtime.rankp1("mpich") )

disp("")
disp('pm.lib.mpi.runtime.rankp1("unknown")')
disp( pm.lib.mpi.runtime.rankp1("unknown") )