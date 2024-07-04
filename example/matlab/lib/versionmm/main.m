cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.lib.versionmm("1.1.1")')
pm.matlab.show( pm.lib.versionmm("1.1.1") )
assert(pm.lib.versionmm("1.1.1") == "1.1.0")

pm.matlab.show()
pm.matlab.show('pm.lib.versionmm("1.1.0")')
pm.matlab.show( pm.lib.versionmm("1.1.0") )
assert(pm.lib.versionmm("1.1.0") == "1.0.0")