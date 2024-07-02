cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.lib.cite()')
disp( pm.lib.cite() )

disp("")
disp('pm.lib.cite([])')
disp( pm.lib.cite([]) )

disp("")
disp('pm.lib.cite("html")')
disp( pm.lib.cite("html") )

disp("")
disp('pm.lib.cite("raw")')
disp( pm.lib.cite("raw") )