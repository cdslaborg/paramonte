cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

cmat = randi([0, 9], 3, 5)
pm.array.verbose(cmat, 1, [0, 2, 1])
pm.array.verbose(cmat, 2, [2, 0, 1, 3, 1])