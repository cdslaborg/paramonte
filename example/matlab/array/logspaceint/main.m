cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

assert(all(pm.array.logspaceint(log(10), log(20)) == [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]));
assert(all(pm.array.logspaceint(log(10), log(20), log(1.5)) == [10, 15]));
assert(all(pm.array.logspaceint(log(10), log(20), [], 2) == [5, 6, 7, 8]));
assert(all(pm.array.logspaceint(log(10), log(100), log(1.5)) == [10, 15, 23, 34, 51, 76]));

disp("pm.array.logspaceint(log(10), log(100), log(1.5))")
disp( pm.array.logspaceint(log(10), log(100), log(1.5)) )