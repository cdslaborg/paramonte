cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.lib.update.version()')
disp( pm.lib.update.version() )

disp("")
disp('pm.lib.update.version([])')
disp( pm.lib.update.version([]) )

disp("")
disp('pm.lib.update.version(false)')
disp( pm.lib.update.version(false) )

disp("")
disp('pm.lib.update.version(true)')
disp( pm.lib.update.version(true) )