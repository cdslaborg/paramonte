cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

assert(pm.array.len(1) == 1);
assert(pm.array.len("") == 0);
assert(pm.array.len([]) == 0);
assert(pm.array.len('paramonte') == 1);
assert(pm.array.len("paramonte") == 1);
assert(pm.array.len(["pm", 'array']) == 2);
assert(pm.array.len(["pm", 'array', []]) == 2);

disp('pm.array.len("paramonte")')
disp( pm.array.len("paramonte") )

disp("pm.array.len('paramonte') == length('paramonte')")
disp( pm.array.len('paramonte') == length('paramonte') )

disp('pm.array.len("paramonte") == length("paramonte")')
disp( pm.array.len("paramonte") == length("paramonte") )