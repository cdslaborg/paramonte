cd(fileparts(mfilename('fullpath'))); % Change working directory to source code directory.
addpath('../../../'); % Add the ParaMonte library root directory to the search path.

pm.matlab.show()
pm.matlab.show('pm.fortran.quote("paramonte")')
pm.matlab.show( pm.fortran.quote("paramonte") )
assert(pm.fortran.quote("paramonte") == """paramonte""")

pm.matlab.show()
pm.matlab.show('pm.fortran.quote(""the" paramonte")')
pm.matlab.show( pm.fortran.quote('"the" paramonte') )
assert(pm.fortran.quote('"the" paramonte') == """""""the"""" paramonte""")