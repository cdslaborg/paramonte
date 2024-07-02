cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.lib.update.weblink()')
disp( pm.lib.update.weblink() )

disp("")
disp('pm.lib.update.weblink([])')
disp( pm.lib.update.weblink([]) )

disp("")
disp('pm.lib.update.weblink(false)')
disp( pm.lib.update.weblink(false) )

disp("")
disp('pm.lib.update.weblink(true)')
disp( pm.lib.update.weblink(true) )