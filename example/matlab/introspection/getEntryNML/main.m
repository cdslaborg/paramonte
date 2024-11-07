cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.introspection.getEntryNML("outputChainSize", 10000, "integer", 1)')
disp( pm.introspection.getEntryNML("outputChainSize", 10000, "integer", 1) )

disp('pm.introspection.getEntryNML("proposalMean", [1, 0], "double", 2)')
disp( pm.introspection.getEntryNML("proposalMean", [1, 0], "double", 2) )

disp('pm.introspection.getEntryNML("mpiEnabled", false, "logical", 1)')
disp( pm.introspection.getEntryNML("mpiEnabled", false, "logical", 1) )

disp('pm.introspection.getEntryNML("description", "This is the description", "string", 1)')
disp( pm.introspection.getEntryNML("description", "This is the description", "string", 1) )

try
    pm.introspection.getEntryNML("varname", "value", "cell", 1);
catch me
    warning(string(me.identifier) + newline + string(me.message));
end

try
    pm.introspection.getEntryNML("varname", [1, 2], "array", 1);
catch me
    warning(string(me.identifier) + newline + string(me.message));
end