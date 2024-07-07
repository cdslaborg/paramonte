cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.stats.dist.cov.getRand(1)')
pm.matlab.show( pm.stats.dist.cov.getRand(1) )

pm.matlab.show()
pm.matlab.show('pm.stats.dist.cov.getRand(3)')
pm.matlab.show( pm.stats.dist.cov.getRand(3) )

pm.matlab.show()
pm.matlab.show('pm.stats.dist.cov.getRand(3, 5)')
pm.matlab.show( pm.stats.dist.cov.getRand(3, 5) )

pm.matlab.show()
pm.matlab.show('pm.stats.dist.cov.getRand(3, [1, 3, 5])')
pm.matlab.show( pm.stats.dist.cov.getRand(3, [1, 3, 5]) )