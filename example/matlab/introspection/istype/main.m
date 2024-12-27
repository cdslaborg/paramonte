cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.introspection.istype(10000, "integer")')
disp( pm.introspection.istype(10000, "integer") )

disp('pm.introspection.istype([1, 0], "double")')
disp( pm.introspection.istype([1, 0], "double") )

disp('pm.introspection.istype(false, "logical")')
disp( pm.introspection.istype(false, "logical") )

disp('pm.introspection.istype("This is the description", "string")')
disp( pm.introspection.istype("This is the description", "string") )

disp('pm.introspection.istype("value", "cell")')
disp( pm.introspection.istype("value", "cell") )

disp('pm.introspection.istype([1, 2], "complex")')
disp( pm.introspection.istype([1, 2], "complex") )

disp('pm.introspection.istype(1 + 2i, "complex")')
disp( pm.introspection.istype(1 + 2i, "complex") )