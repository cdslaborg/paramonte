cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.str.index("paramonte", "")')
                pm.str.index("paramonte", "")
assert(pm.str.index("paramonte", "") == 0)

pm.matlab.show()
pm.matlab.show('pm.str.index("paramonte", "M")')
                pm.str.index("paramonte", "M")
assert(pm.str.index("paramonte", "M") == 0)

pm.matlab.show()
pm.matlab.show('pm.str.index("paramonte", "monte")')
                pm.str.index("paramonte", "monte")
assert(pm.str.index("paramonte", "monte") == 5)