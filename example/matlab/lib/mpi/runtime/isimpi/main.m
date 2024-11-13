cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

disp('itis = pm.lib.mpi.runtime.isimpi()')
      itis = pm.lib.mpi.runtime.isimpi()
disp("")