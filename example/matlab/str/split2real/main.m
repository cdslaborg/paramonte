cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.str.split2real("3.1.2")')
pm.matlab.show( pm.str.split2real("3.1.2") )
assert(all(pm.str.split2real("3.1.2") == [3, 1, 2]))

pm.matlab.show()
pm.matlab.show('pm.str.split2real("3-1-2.5", "-")')
pm.matlab.show( pm.str.split2real("3-1-2.5", "-") )
assert(all(pm.str.split2real("3-1-2.5", "-") == [3, 1, 2.5]))