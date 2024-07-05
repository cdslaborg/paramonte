cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.str.isnumeric("paramonte")')
                pm.str.isnumeric("paramonte")
assert(~pm.str.isnumeric("paramonte"))

pm.matlab.show()
pm.matlab.show('pm.str.isnumeric("1")')
                pm.str.isnumeric("1")
assert(pm.str.isnumeric("1"))

pm.matlab.show()
pm.matlab.show('pm.str.isnumeric("1.2")')
                pm.str.isnumeric("1.2")
assert(pm.str.isnumeric("1.2"))

pm.matlab.show()
pm.matlab.show('pm.str.isnumeric("1.2 + 1i")')
                pm.str.isnumeric("1.2 + 1i")
assert(pm.str.isnumeric("1.2 + 1i"))

pm.matlab.show()
pm.matlab.show('pm.str.isnumeric("[1.2, 1i]")')
                pm.str.isnumeric("[1.2, 1i]")
assert(~pm.str.isnumeric("[1.2, 1i]"))