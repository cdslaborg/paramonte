cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show("pm.math.mincirc(5)")
pm.matlab.show( pm.math.mincirc(5) )

pm.matlab.show()
pm.matlab.show("pm.math.mincirc(21)")
pm.matlab.show( pm.math.mincirc(21) )

[nrow, ncol] = pm.math.mincirc(5);
assert(all([nrow, ncol] == [3, 2]))

[nrow, ncol] = pm.math.mincirc(5);
assert(all([nrow, ncol] == [7, 3]))