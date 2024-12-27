cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.introspection.verified(10000, "integer", 1)')
disp( pm.introspection.verified(10000, "integer", 1) )

disp('pm.introspection.verified([10000, 12], "integer", 1)')
disp( pm.introspection.verified([10000, 12], "integer", 1) )

disp('pm.introspection.verified([10000, 12], "integer", 2)')
disp( pm.introspection.verified([10000, 12], "integer", 2) )

disp('pm.introspection.verified([1, 0], "double", 1)')
disp( pm.introspection.verified([1, 0], "double", 1) )

disp('pm.introspection.verified([1, 0], "double", 2)')
disp( pm.introspection.verified([1, 0], "double", 2) )

disp('pm.introspection.verified([1, 0], "double", 3)')
disp( pm.introspection.verified([1, 0], "double", 3) )

disp('pm.introspection.verified(false, "logical", 2)')
disp( pm.introspection.verified(false, "logical", 2) )

disp('pm.introspection.verified("This is the description", "string", 0)')
disp( pm.introspection.verified("This is the description", "string", 0) )

disp('pm.introspection.verified("This is the description", "string", 1)')
disp( pm.introspection.verified("This is the description", "string", 1) )

disp('pm.introspection.verified("value", "cell", 1)')
disp( pm.introspection.verified("value", "cell", 1) )

disp('pm.introspection.verified([1, 2], "complex", 1)')
disp( pm.introspection.verified([1, 2], "complex", 1) )

disp('pm.introspection.verified([1, 2], "complex", 2)')
disp( pm.introspection.verified([1, 2], "complex", 2) )

disp('pm.introspection.verified(1 + 2i, "complex", 1)')
disp( pm.introspection.verified(1 + 2i, "complex", 1) )