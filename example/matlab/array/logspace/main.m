cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.array.logspace(log(10), log(20))
pm.array.logspace(log(10), log(20), log(1.5)) % = 10.000000000000002  15.000000000000007
pm.array.logspace(log(10), log(20), [], 2)
pm.array.logspace(log(10), log(20), log(1.5), 3) % = 12.549091579662001  19.591528026787028