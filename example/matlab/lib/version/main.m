cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp("")
disp('pm.lib.version()')
disp( pm.lib.version() )

disp("")
disp('pm.lib.version("1.1.1")') % 1.1.0
disp( pm.lib.version("1.1.1") ) % 1.1.0

disp("")
disp('pm.lib.version("1.1.0")') % 1.0.0
disp( pm.lib.version("1.1.0") ) % 1.0.0

assert(pm.lib.version("1.1.1") == "1.1.0")
assert(pm.lib.version("1.1.0") == "1.0.0")