cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

disp('[nproc, rankp1] = pm.lib.mpi.runtime.hydra()')
      [nproc, rankp1] = pm.lib.mpi.runtime.hydra()
disp("")