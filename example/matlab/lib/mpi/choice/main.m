cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.lib.mpi.choice()')
disp( pm.lib.mpi.choice() )