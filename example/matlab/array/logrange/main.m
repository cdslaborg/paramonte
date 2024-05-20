cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

assert(all(pm.array.logrange(1, 1000, 10) == [1, 2, 4, 8, 16, 32, 63, 126, 251, 501, 1000]));
assert(all(pm.array.logrange(1, 1000, 20) == [1, 2, 3, 4, 6, 8, 11, 16, 22, 32, 45, 63, 89, 126, 178, 251, 355, 501, 708, 1000]));
pm.array.logrange(1, 1000, [])
pm.array.logrange(1, 100)