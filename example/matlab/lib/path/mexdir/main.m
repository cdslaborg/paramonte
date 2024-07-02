cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.lib.path.mexdir()')
disp( pm.lib.path.mexdir() )

disp("")
disp('pm.lib.path.mexdir("sampling")')
disp( pm.lib.path.mexdir("sampling") )

disp("")
disp('pm.lib.path.mexdir("sampling", "release")')
disp( pm.lib.path.mexdir("sampling", "release") )

disp("")
disp('pm.lib.path.mexdir("sampling", ["release", "GNU"])')
disp( pm.lib.path.mexdir("sampling", ["release", "GNU"]) )

disp("")
disp('pm.lib.path.mexdir("sampling", ["native", "gnu"])')
disp( pm.lib.path.mexdir("sampling", ["native", "gnu"]) )

disp("")
disp('pm.lib.path.mexdir("sampling", ["native", "GNU"])')
disp( pm.lib.path.mexdir("sampling", ["native", "GNU"]) )

disp("")
disp('pm.lib.path.mexdir("sampling", ["native", "Intel"])')
disp( pm.lib.path.mexdir("sampling", ["native", "Intel"]) )
