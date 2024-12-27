cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.introspection.islenleq(10000, 1)')
disp( pm.introspection.islenleq(10000, 1) )

disp('pm.introspection.islenleq([1, 0], 2)')
disp( pm.introspection.islenleq([1, 0], 2) )

disp('pm.introspection.islenleq(false, 1)')
disp( pm.introspection.islenleq(false, 1) )

disp('pm.introspection.islenleq("This is the description", 1)')
disp( pm.introspection.islenleq("This is the description", 1) )

disp('pm.introspection.islenleq("value", 1)')
disp( pm.introspection.islenleq("value", 1) )

disp('pm.introspection.islenleq([1, 2], 1)')
disp( pm.introspection.islenleq([1, 2], 1) )

disp('pm.introspection.islenleq(1 + 2i, 1)')
disp( pm.introspection.islenleq(1 + 2i, 1) )