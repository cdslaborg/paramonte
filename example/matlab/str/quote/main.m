cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.str.quote("paramonte")')
pm.matlab.show( pm.str.quote("paramonte") )
assert(pm.str.quote("paramonte") == """paramonte""")

pm.matlab.show()
pm.matlab.show('pm.str.quote(""the" paramonte")')
pm.matlab.show( pm.str.quote('"the" paramonte') )
assert(pm.str.quote('"the" paramonte') == """""""the"""" paramonte""")