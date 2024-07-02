cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.lib.path.mex()')
disp( pm.lib.path.mex() )

disp("")
disp('pm.lib.path.mex("sampling")')
disp( pm.lib.path.mex("sampling") )

disp("")
disp('pm.lib.path.mex("sampling", "release")')
disp( pm.lib.path.mex("sampling", "release") )

disp("")
disp('pm.lib.path.mex("sampling", ["release", "gnu"])')
disp( pm.lib.path.mex("sampling", ["release", "gnu"]) )

disp("")
disp('pm.lib.path.mex("sampling", ["release", "GNU"])')
disp( pm.lib.path.mex("sampling", ["release", "GNU"]) )

disp("")
disp('pm.lib.path.mex("sampling", ["native", "gnu"])')
disp( pm.lib.path.mex("sampling", ["native", "gnu"]) )

disp("")
disp('pm.lib.path.mex("sampling", ["native", "GNU"])')
disp( pm.lib.path.mex("sampling", ["native", "GNU"]) )

disp("")
disp('pm.lib.path.mex("sampling", ["native", "Intel"])')
disp( pm.lib.path.mex("sampling", ["native", "Intel"]) )
