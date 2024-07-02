cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.introspection.istype(10000, "integer", 1)')
disp( pm.introspection.istype(10000, "integer", 1) )

disp('pm.introspection.istype([1, 0], "double", 2)')
disp( pm.introspection.istype([1, 0], "double", 2) )

disp('pm.introspection.istype(false, "logical", 1)')
disp( pm.introspection.istype(false, "logical", 1) )

disp('pm.introspection.istype("This is the description", "string", 1)')
disp( pm.introspection.istype("This is the description", "string", 1) )

disp('pm.introspection.istype("value", "cell", 1)')
disp( pm.introspection.istype("value", "cell", 1) )

disp('pm.introspection.istype([1, 2], "complex", 1)')
disp( pm.introspection.istype([1, 2], "complex", 1) )

disp('pm.introspection.istype(1 + 2i, "complex", 1)')
disp( pm.introspection.istype(1 + 2i, "complex", 1) )