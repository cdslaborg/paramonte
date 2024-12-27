cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

disp('pm.fortran.getEntryNML("outputChainSize", 10000, "integer", 1)')
disp( pm.fortran.getEntryNML("outputChainSize", 10000, "integer", 1) )

disp('pm.fortran.getEntryNML("proposalMean", [1, 0], "double", 2)')
disp( pm.fortran.getEntryNML("proposalMean", [1, 0], "double", 2) )

disp('pm.fortran.getEntryNML("mpiEnabled", false, "logical", 1)')
disp( pm.fortran.getEntryNML("mpiEnabled", false, "logical", 1) )

disp('pm.fortran.getEntryNML("description", "This is the description", "string", 1)')
disp( pm.fortran.getEntryNML("description", "This is the description", "string", 1) )

try
    pm.fortran.getEntryNML("varname", "value", "cell", 1);
catch me
    warning(string(me.identifier) + newline + string(me.message));
end

try
    pm.fortran.getEntryNML("varname", [1, 2], "array", 1);
catch me
    warning(string(me.identifier) + newline + string(me.message));
end